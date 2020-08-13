include("startup_loadpackages.jl")
include("accessories.jl")

config = parsefile("Data.toml")
allmeta = CSV.File(config["tables"]["joined_metadata"], pool=false) |> DataFrame
kidsmeta = filter(row-> !startswith(row.sample, "M") && !ismissing(row.correctedAgeDays), allmeta)

## Feature tables
@warn "Loading feature tables"

species = widen2comm(taxonomic_profiles(filefilter=f-> sampleid(stoolsample(basename(f))) in kidsmeta.sample)...)
# Total sum scaling - function in Microbiome
relativeabundance!(species)

kos = widen2comm(functional_profiles(kind="kos_names_relab", filefilter=f-> sampleid(stoolsample(basename(f))) in kidsmeta.sample)..., featurecol=:func)
kos = view(kos, species=map(x-> !in(x, ("UNMAPPED", "UNGROUPED")), featurenames(kos)))
stratkos = widen2comm(functional_profiles(kind="kos_names_relab", filefilter=f-> sampleid(stoolsample(basename(f))) in kidsmeta.sample, stratified=true)..., featurecol=:func)
stratkos = view(stratkos, species=map(x-> !in(x, ("UNMAPPED", "UNGROUPED")), featurenames(stratkos)))
pfams = widen2comm(functional_profiles(kind="pfams_names_relab", filefilter=f-> sampleid(stoolsample(basename(f))) in kidsmeta.sample)..., featurecol=:func)
pfams = view(pfams, species=map(x-> !in(x, ("UNMAPPED", "UNGROUPED")), featurenames(pfams)))
ecs = widen2comm(functional_profiles(kind="ecs_names_relab", filefilter=f-> sampleid(stoolsample(basename(f))) in kidsmeta.sample)..., featurecol=:func)
ecs = view(ecs, species=map(x-> !in(x, ("UNMAPPED", "UNGROUPED")), featurenames(ecs)))

unirefs = widen2comm(functional_profiles(kind="genefamilies_relab", filefilter=f-> sampleid(stoolsample(basename(f))) in kidsmeta.sample)..., featurecol=:func)

@warn "Getting metadata and subgroups"

### Just get metadata found in tax/func profiles, and in same order

allmeta.ageLabel = map(eachrow(allmeta)) do row
    ismissing(row.correctedAgeDays) && return missing
    row.correctedAgeDays < 365 && return "1 and under"
    row.correctedAgeDays < 365*2 && return "1 to 2"
    return "2 and over"
end
dropmissing!(allmeta, :ageLabel)

filter!(allmeta) do row
    !ismissing(row.cogScore) ||
    !ismissing(row.hires_total)
end

# make sure all tables have the same samples in the same order
allsamples = intersect(allmeta.sample, map(sitenames, (species, unirefs, ecs, kos, pfams))...) |> collect |> sort
filter!(row-> row.sample in allsamples, allmeta)

species  = view(species, sites=allmeta.sample)  |> copy
unirefs  = view(unirefs, sites=allmeta.sample)  |> copy
kos      = view(kos, sites=allmeta.sample)      |> copy
stratkos = view(stratkos, sites=allmeta.sample) |> copy
pfams    = view(pfams, sites=allmeta.sample)    |> copy
ecs      = view(ecs, sites=allmeta.sample)      |> copy

@assert samplenames(species) == samplenames(unirefs) == samplenames(pfams) == samplenames(kos) == samplenames(ecs)

## Subgroup indexes

ukids, oldkids = let samples = Set(sampleid.(uniquetimepoints(allmeta.sample, takefirst=true, samplefilter=iskid)))
    (map(s-> in(s, samples), allmeta.sample),
    map(row-> in(row.sample, samples) && row.ageLabel != "1 and under", eachrow(allmeta)))
end

@warn "Getting accessory genes"

Microbiome.prevalence(a, minabundance::Float64=0.0001) = mean(x-> present(x, minabundance), (y for y in a))

unirefprevfilt = let c = 0
    filt = map(eachrow(occurrences(unirefs))) do row
        c % 10_000 == 0 && @info "    processed $c genes"
        u1_prev = prevalence(row[map(x->  x == "1 and under", allmeta.ageLabel)], 0.)
        o1_prev = prevalence(row[map(x->  x != "1 and under", allmeta.ageLabel)], 0.)
        c+=1
        (any(>(0.05), [u1_prev, o1_prev]), all(<(0.9), [u1_prev, o1_prev]) )
    end
end
unirefprevalent = view(unirefs, species=[p[1] for p in unirefprevfilt])
unirefaccessory = view(unirefs, species=[p[2] for p in unirefprevfilt])

### Additional info
@warn "Adding additional info to metadata tables"

### Calculate some additional metadata

allmeta.shannon = shannon(species)

unmappedidx = findall(u-> occursin("UNMAPPED",u), featurenames(unirefs))[1]

allmeta.identifiable_unirefs = 1 .- Vector(occurrences(unirefs)[unmappedidx,:])

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

ukidsmeta = view(allmeta, ukids, :)
oldkidsmeta = view(allmeta, oldkids, :)

@warn "Done!"
