# NOTE: This can take 30 minutes to a few hours to run

@warn "Loading packages"

ENV["GKSwstype"] = "100"
using ECHOAnalysis
using DataFrames
using SQLite
using StatsPlots
using PrettyTables
using CSV
using Microbiome
using Distances
using MultivariateStats
using Pkg.TOML: parsefile
using Clustering
using Combinatorics
using MicrobiomePlots
using BiobakeryUtils
using MultipleTesting
using ProgressMeter
using JLD2


rounder = Dict(0 => (v,i) -> typeof(v) <: AbstractFloat ? round(v,digits=3) : v)
# print ~15 random rows
randrowfilter(data, i) = rand() < (1 / size(data, 1)) * 15
@ptconfclean # clear any previous configuration
@ptconf formatter = rounder nosubheader=true screen_size=(20,120) filters_row=(randrowfilter,)

include("accessories.jl")

config = parsefile("data/data.toml")
allmeta = getmgxmetadata()

allmeta = getmgxmetadata(samples=uniquetimepoints(allmeta.sample, takefirst=false))

allmeta.breastfeeding = map(eachrow(allmeta)) do row
    ismissing(row.breastFedPercent) && return missing
    !(row.breastFedPercent isa Number) && error(":breastFedPercent should be a number or missing")
    if row.breastFedPercent < 5
        return "exclussive formula"
    elseif row.breastFedPercent > 80
        return "exclussive breast"
    else
        return "mixed"
    end
end

allmeta.childWeight = map(w-> ismissing(w) ? w : parse(Float64, w), allmeta.childWeight)
allmeta.BMI_calc = map(row-> row.childWeight / (row.childHeight^2) * 703, eachrow(allmeta))

figures = config["output"]["figures"]
figures = joinpath(figures, "figure1")
tables = config["output"]["tables"]
isdir(figures) || mkpath(figures)
isdir(tables) || mkpath(tables)

## Feature tables
@warn "Loading feature tables"

### Load taxonomic profiles from SQLite database
taxdb = SQLite.DB(config["sqlite"]["taxa"]["path"])
### Function in ECHOAnalysis package
species = sqlprofile(taxdb, tablename="taxa", kind="species")
### Keep only samples in the metadata
species = view(species, sites=allmeta.sample) |> copy
### Total sum scaling - function in Microbiome
relativeabundance!(species)

### Same steps for functional profiles
unirefdb = SQLite.DB(config["sqlite"]["uniref90"]["path"])
unirefs = sqlprofile(unirefdb, tablename="genefamilies_relab", kind="genefamilies_relab")
### We have a couple of taxonomic profiles that don't have functional profiles,
### this gets tax profiles in the same order
species = view(species, sites=sitenames(unirefs)) |> copy

kodb = SQLite.DB(config["sqlite"]["ko"]["path"])
kos = sqlprofile(kodb, tablename="ko_names_relab", kind="ko_names_relab")
kos = view(kos, species=map(x-> !in(x, ("UNMAPPED", "UNGROUPED")), featurenames(kos)))
pfamdb = SQLite.DB(config["sqlite"]["pfam"]["path"])
pfams = sqlprofile(pfamdb, tablename="pfam_names_relab", kind="pfam_names_relab")
pfams = view(pfams, species=map(x-> !in(x, ("UNMAPPED", "UNGROUPED")), featurenames(pfams)))
ecdb = SQLite.DB(config["sqlite"]["ec"]["path"])
ecs = sqlprofile(ecdb, tablename="ec_names_relab", kind="ec_names_relab")
ecs = view(ecs, species=map(x-> !in(x, ("UNMAPPED", "UNGROUPED")), featurenames(ecs)))

species = view(species, sites=sitenames(ecs)) |> copy

@warn "Getting metadata and subgroups"

### double check that everything has samples in the same order for later indexing
### `@assert` will throw an error if expression is not `true`
@assert samplenames(species) == samplenames(unirefs) == samplenames(pfams) == samplenames(kos) == samplenames(ecs)
### Just get metadata found in tax/func profiles, and in same order
allmeta = getmgxmetadata(samples=sitenames(species))
dropmissing!(allmeta, :ageLabel)
species = view(species, sites=allmeta.sample) |> copy
unirefs = view(unirefs, sites=allmeta.sample) |> copy
kos     = view(kos, sites=allmeta.sample)     |> copy
pfams   = view(pfams, sites=allmeta.sample)   |> copy
ecs     = view(ecs, sites=allmeta.sample)     |> copy

## Subgroup indexes

allmoms = let samples = Set(sampleid.(uniquetimepoints(allmeta.sample, takefirst=false, samplefilter=ismom)))
    map(s-> in(s, samples), allmeta.sample)
end

allkids = let samples = Set(sampleid.(uniquetimepoints(allmeta.sample, takefirst=false, samplefilter=iskid)))
    map(s-> in(s, samples), allmeta.sample)
end

umoms = let samples = Set(sampleid.(uniquetimepoints(allmeta.sample, takefirst=true, samplefilter=ismom)))
    map(s-> in(s, samples), allmeta.sample)
end

ukids, oldkids = let samples = Set(sampleid.(uniquetimepoints(allmeta.sample, takefirst=true, samplefilter=iskid)))
    (map(s-> in(s, samples), allmeta.sample),
    map(row-> in(row.sample, samples) && row.ageLabel == "2 and over", eachrow(allmeta)))
end

uboth = map(any, zip(umoms, ukids))

### Metadata subgroups

allmomsmeta = view(allmeta, allmoms, :)
umomsmeta = view(allmeta, umoms, :)
allkidsmeta = view(allmeta, allkids, :)
ukidsmeta = view(allmeta, ukids, :)
oldkidsmeta = view(allmeta, oldkids, :)
ubothmeta = view(allmeta, uboth, :)

@warn "Getting accessory genes"

Microbiome.prevalence(a, minabundance::Float64=0.0001) = mean(x-> present(x, minabundance), (y for y in a))

unirefprevfilt = let c = 0
    filt = map(eachrow(occurrences(unirefs))) do row
        c % 10_000 == 0 && @info "    processed $c genes"
        u1_prev = prevalence(row[map(x->  x == "1 and under", allmeta.ageLabel)], 0.)
        o1_prev = prevalence(row[map(x->  x != "mom" && x != "1 and under", allmeta.ageLabel)], 0.)
        mom_prev = prevalence(row[map(x-> x == "mom", allmeta.ageLabel)], 0.)
        c+=1
        (any(>(0.05), [u1_prev, o1_prev, mom_prev]), all(<(0.9), [u1_prev, o1_prev, mom_prev]) )
    end
end
unirefprevalent = view(unirefs, species=[p[1] for p in unirefprevfilt])
unirefaccessory = view(unirefs, species=[p[2] for p in unirefprevfilt])

### Additional info
@warn "Adding additional info to metadata tables"

allmeta.shannon = shannon(species)
allmeta.identifiable_unirefs = 1 .- Vector(occurrences(unirefs)[1,:])
allmeta.n_unirefs = vec(sum(!=(0.), occurrences(unirefs)[1:end,:], dims=1))
allmeta.subject_type = [ismom(s) ? "Mother" :
                        iskid(s) ? "Child"  :
                        error("Not mom or kid: $s") for s in allmeta.sample]

allmeta.lowres_total = map(row-> sum(row[[:white_matter_volume, :grey_matter_volume, :csf_volume]]), eachrow(allmeta))
allmeta.white_matter_normed = allmeta.white_matter_volume ./ allmeta.lowres_total
allmeta.grey_matter_normed = allmeta.grey_matter_volume ./ allmeta.lowres_total
allmeta.csf_normed = allmeta.csf_volume ./ allmeta.lowres_total

allmeta.cerebellar_normed = allmeta.cerebellar ./ allmeta.hires_total
allmeta.neocortical_normed = allmeta.neocortical ./ allmeta.hires_total
allmeta.subcortical_normed = allmeta.subcortical ./ allmeta.hires_total
allmeta.limbic_normed = allmeta.limbic ./ allmeta.hires_total

@warn "Done!"
