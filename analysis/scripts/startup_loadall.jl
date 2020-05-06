# NOTE: This can take 30 minutes to a few hours to run

include("startup_loadpackages.jl")
include("accessories.jl")

config = parsefile("Data.toml")
allmeta = CSV.File(config["tables"]["joined_metadata"]) |> DataFrame

## Feature tables
@warn "Loading feature tables"

species = widen2comm(taxonomic_profiles()...)
# Total sum scaling - function in Microbiome
relativeabundance!(species)

# all but unirefs
include("startup_loadfunctional.jl")
unirefs = widen2comm(functional_profiles(kind="genefamilies_relab")..., featurecol=:func)

@warn "Getting metadata and subgroups"

### Just get metadata found in tax/func profiles, and in same order

allsamples = intersect(map(sitenames, (species, unirefs, ecs, kos, pfams))...) |> collect |> sort

allmeta.ageLabel = map(eachrow(allmeta)) do row
    startswith(row.sample, "M") && return "mom"
    ismissing(row.correctedAgeDays) && return missing
    row.correctedAgeDays < 365 && return "1 and under"
    row.correctedAgeDays < 365*2 && return "1 to 2"
    return "2 and over"
end

dropmissing!(allmeta, :ageLabel)

species = view(species, sites=allmeta.sample) |> copy
unirefs = view(unirefs, sites=allmeta.sample) |> copy
kos     = view(kos, sites=allmeta.sample)     |> copy
pfams   = view(pfams, sites=allmeta.sample)   |> copy
ecs     = view(ecs, sites=allmeta.sample)     |> copy

@assert samplenames(species) == samplenames(unirefs) == samplenames(pfams) == samplenames(kos) == samplenames(ecs)

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

### Calculate some additional metadata

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

### Metadata subgroups

allmomsmeta = view(allmeta, allmoms, :)
umomsmeta = view(allmeta, umoms, :)
allkidsmeta = view(allmeta, allkids, :)
ukidsmeta = view(allmeta, ukids, :)
oldkidsmeta = view(allmeta, oldkids, :)
ubothmeta = view(allmeta, uboth, :)

@warn "Done!"
