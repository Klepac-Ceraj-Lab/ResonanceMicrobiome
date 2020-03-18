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
using Colors

kodb = SQLite.DB("~/Desktop/ko_profiles.sqlite")
allmeta = getmgxmetadata("/Users/ksb/Desktop/metadata.sqlite")
allmeta = getmgxmetadata("/Users/ksb/Desktop/metadata.sqlite",samples=uniquetimepoints(allmeta.sample, takefirst=false))
include("accessories.jl")
# add_functional_profiles(kodb, "data/engaging", stratified=true, kind="kos_relab")

koslong = DataFrame(DBInterface.execute(kodb, "SELECT * FROM kos_relab"))
nko = get_neuroactive_kos()

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
koslong.agelabel = [agelabels[s] for s in koslong.sample]

(s,t) = let geneset = nko["Glutamate degradation"]
    filt = filter(row-> row.ko in geneset, koslong)
    samples = by(filt, [:sample, :taxon]) do sample
        (total=sum(sample.abundance), agelabel=first(sample).agelabel)
    end

    sort!(samples, [:agelabel, :sample])
    totals = filter(row-> isnothing(row.taxon), samples)
    filter!(row-> !isnothing(row.taxon), samples)
    totals.ord = invperm(sortperm(totals, [:agelabel, :total]))
    alltaxa = Set(samples.taxon)

    for s in unique(samples.sample)
        missing_tax = setdiff(alltaxa, samples[samples.sample .== s, :taxon])
        for t in missing_tax
            push!(samples, (sample=s, taxon=t, total=0., agelabel=agelabels[s]))
        end
    end
    xs = dictionary(s=> i for (s,i) in eachrow(totals[!,[:sample,:ord]]))
    samples.x = [xs[s] for s in samples.sample]
    sort!(samples, :x)
    (samples, totals)
end

CSV.write("/Users/ksb/Desktop/samples.csv", s)

##

scene, layout = layoutscene()
ax = layout[1,1] = LAxis(scene)
p1 = barplot!(ax, Position.stack, Data(s), Group(:taxon), :x, :total)

tightlimits!(ax)
hidexdecorations!(ax)

annbar = layout[2,1] = LAxis(scene, height=100)
linkxaxes!(ax, annbar)
tightlimits!(annbar)
hidexdecorations!(annbar)
hideydecorations!(annbar)

for (l, c) in zip(unique(s.agelabel), ColorBrewer.palette("Set3", 4))
    (start, stop) = extrema(s[s.agelabel .== l, :x])
    poly!(annbar, Point2f0[(start-0.5,0), (stop+0.5, 0), (stop+0.5, 1), (start-0.5, 1)], color = c)
end

p1legend= layout[1,2] = LLegend(scene, height=Auto(true), width=Auto(false),
                            titlevisible=false, patchcolor=:transparent)

for (t,p) in zip(unique(s.taxon), p1.plots)
    push!(p1legend,
        LegendEntry(t, MarkerElement(color = p.attributes.color[], marker = :rect, strokecolor = :black))
        )
end

scene


### incomplete

# function plotstrat(genedf, geneset; samplecol=:sample, groupcol=:taxon,
#                                     anncols = [:agelabel], abundcol=:abundance, groupn=5)
#
#     annotations = dictionary(a => )
#
#     filt = filter(row-> row.ko in geneset, genedf)
#     samples = by(filt, [samplecol, groupcol]) do sample
#             DataFrame(
#                 :total=>sum(sample[!, abundcol])
#                 [a => first(sample)[a] for a in anncols]...)
#     end
#
#     sort!(samples, [anncols..., samplecol])
#     totals = filter(row-> isnothing(row[groupcol]), samples)
#     filter!(row-> !isnothing(row[groupcol]), samples)
#     totals.ord = invperm(sortperm(totals, [anncols..., :total]))
#     alltaxa = Set(samples.taxon)
#
#     for s in unique(samples[!,samplecol])
#         missing_groups = setdiff(alltaxa, samples[samples[!,samplecol] .== g, groupcol])
#         for mg in missing_groups
#             push!(samples, (; samplcol=s, groupcol=mg, total=0., agelabel=agelabels[s]))
#         end
#     end
#     xs = dictionary(s=> i for (s,i) in eachrow(totals[!,[:sample,:ord]]))
#     samples.x = [xs[s] for s in samples.sample]
#     sort!(samples, :x)
#     (samples, totals)
# end
