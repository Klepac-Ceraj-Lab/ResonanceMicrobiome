using Pkg
Pkg.activate(".")
using ECHOAnalysis
using DataFrames
using CSV
using SQLite
using Microbiome
using PrettyTables
import Pkg.TOML: parsefile
using MakieLayout
using Makie
using StatsMakie

config = parsefile("data/data.toml")
kodb = SQLite.DB(config["sqlite"]["ko"]["path"])
include("startup_loadmetadata.jl")
include("accessories.jl")
add_functional_profiles(kodb, "data/engaging", stratified=true, kind="kos_relab")

koslong = DataFrame(DBInterface.execute(kodb, "SELECT * FROM kos_relab"))
nko = get_neuroactive_kos()

ko1 = nko["Kynurenine degradation"][1]

# replace this with `select!()` once it's added https://github.com/JuliaData/DataFrames.jl/pull/2080
koslong = hcat(koslong, DataFrame(map(eachrow(koslong)) do row
    m = match(r"^(\w+)\|?(?:(?:g__([\w]+)\.s__([\w]+))?|(unclassified))$", row.function)
    (ko, genus, species, unclass) = m.captures
    if isnothing(genus) && isnothing(unclass)
        ko == row.function || error("Weird function $(row.function)")
        taxon = nothing
    else
        isnothing(unclass) ? taxon = "$genus $species" : taxon = "unclassified"
    end

    (ko=ko, taxon=taxon)
    end
))
filter!(row-> !ismissing(row.ageLabel), allmeta)
let samples = Set(allmeta.sample)
    filter!(row-> row.sample in samples, koslong)
end
agelabels = dictionary(sample=>age for (sample, age) in eachrow(allmeta[!, [:sample, :ageLabel]]))

let geneset = nko["Glutamate degradation"]
    filt = filter(row-> row.ko in geneset, koslong)
    samples = by(filt, [:sample, :taxon]) do sample
        (;total=sum(sample.abundance))
    end
    samples.ageLabel = [agelabels[s] for s in samples.sample]
    sort!(samples, [:ageLabel, :sample])
    totals = filter(row-> isnothing(row.taxon), samples)
    totals.ord = sortperm(totals.total)
    x = dictionary(s=> i for (s,i) in eachrow(totals[!,[:sample,:ord]]))
end

