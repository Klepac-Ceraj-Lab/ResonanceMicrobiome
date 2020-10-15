include("startup_loadpackages.jl")

config = parsefile("Data.toml")
allmeta = CSV.File(config["tables"]["joined_metadata"], pool=false) |> DataFrame

filter!(:correctedAgeDays=> !ismissing, allmeta)

let samples = Set(sampleid.(uniquetimepoints(stoolsample.(allmeta.sample), takefirst=false, samplefilter=iskid)))
    filter!(:sample => s-> s ∈ samples, allmeta)  
end

allmeta.simple_race = map(allmeta.simple_race) do r
    (ismissing(r) || occursin("Decline", r) || occursin("Unknown", r)) && return missing
    occursin("Mixed", r) && return "Mixed"
    occursin('\n', r) && error(":simple_race has newlines")
    return r
end

## Sanity checks

@assert setdiff(allmeta.simple_race, 
                Set(["Caucasian / White", 
                     "Mixed", 
                     "Asian ",
                     missing,
                     "African American / Black", 
                     "Native American / Alaskan Native"])
                ) |> length == 0


## Feature tables
@warn "Loading feature tables"
species = widen2comm(taxonomic_profiles(filefilter=f-> sampleid(stoolsample(basename(f))) in allmeta.sample)...)

missingsp = setdiff(allmeta.sample, samplenames(species))

# Total sum scaling - function in Microbiome
relativeabundance!(species)

kos = widen2comm(functional_profiles(kind="kos_names_relab", filefilter=f-> sampleid(stoolsample(basename(f))) in allmeta.sample)..., featurecol=:func)
kos = view(kos, species=map(x-> !in(x, ("UNMAPPED", "UNGROUPED")), featurenames(kos))) |> copy
stratkos = widen2comm(functional_profiles(kind="kos_names_relab", filefilter=f-> sampleid(stoolsample(basename(f))) in allmeta.sample, stratified=true)..., featurecol=:func)
stratkos = view(stratkos, species=map(x-> !in(x, ("UNMAPPED", "UNGROUPED")), featurenames(stratkos))) |> copy
pfams = widen2comm(functional_profiles(kind="pfams_names_relab", filefilter=f-> sampleid(stoolsample(basename(f))) in allmeta.sample)..., featurecol=:func)
pfams = view(pfams, species=map(x-> !in(x, ("UNMAPPED", "UNGROUPED")), featurenames(pfams))) |> copy
ecs = widen2comm(functional_profiles(kind="ecs_names_relab", filefilter=f-> sampleid(stoolsample(basename(f))) in allmeta.sample)..., featurecol=:func)
ecs = view(ecs, species=map(x-> !in(x, ("UNMAPPED", "UNGROUPED")), featurenames(ecs))) |> copy

unirefs = widen2comm(functional_profiles(kind="genefamilies_relab", filefilter=f-> sampleid(stoolsample(basename(f))) in allmeta.sample)..., featurecol=:func)


missingko = setdiff(allmeta.sample, samplenames(kos))
missingpfam = setdiff(allmeta.sample, samplenames(pfams))
missingec = setdiff(allmeta.sample, samplenames(ecs))
missingur = setdiff(allmeta.sample, samplenames(unirefs))
@assert missingko == missingpfam == missingec == missingur
@assert length(setdiff(missingsp, missingko)) == 0
open("/lovelace/echo/analysis/bb3_missing.txt", "w") do io
    println.(Ref(io), missingko)
end

@warn "Getting metadata and subgroups"

### Just get metadata found in tax/func profiles, and in same order

allmeta.ageLabel = map(eachrow(allmeta)) do row
    ismissing(row.correctedAgeDays) && error("No age for $(row.sample)")
    row.correctedAgeDays < 365 && return "1 and under"
    row.correctedAgeDays < 365*2 && return "1 to 2"
    return "2 and over"
end

filter!([:cogScore, :braintotal]=> (cs, bt)-> any(!ismissing, (cs, bt)), allmeta)

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

unirefprevalent, unirefaccessory = accessorygenes(unirefs, [
    row-> prevalence(row[map(x->  x == "1 and under", allmeta.ageLabel)], 0.),
    row-> prevalence(row[map(x->  x != "1 and under", allmeta.ageLabel)], 0.)
    ],
    lower=0.1, upper=0.9)

pfamprevalent, pfamaccessory = accessorygenes(pfams, [
    row-> prevalence(row[map(x->  x == "1 and under", allmeta.ageLabel)], 0.),
    row-> prevalence(row[map(x->  x != "1 and under", allmeta.ageLabel)], 0.)
    ],
    lower=0.1, upper=0.9
)

koprevalent, koaccessory = accessorygenes(kos, [
    row-> prevalence(row[map(x->  x == "1 and under", allmeta.ageLabel)], 0.),
    row-> prevalence(row[map(x->  x != "1 and under", allmeta.ageLabel)], 0.)
    ],
    lower=0.1, upper=0.9
)

ecprevalent, ecaccessory = accessorygenes(ecs, [
    row-> prevalence(row[map(x->  x == "1 and under", allmeta.ageLabel)], 0.),
    row-> prevalence(row[map(x->  x != "1 and under", allmeta.ageLabel)], 0.)
    ],
    lower=0.1, upper=0.9
)


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


allmeta.white_matter_normed    = allmeta.white_matter    ./ allmeta.braintotal
allmeta.gray_matter_normed     = allmeta.gray_matter     ./ allmeta.braintotal
allmeta.csf_normed             = allmeta.csf             ./ allmeta.braintotal
allmeta.hippocampus_normed     = allmeta.hippocampus     ./ allmeta.braintotal
allmeta.thalamus_normed        = allmeta.thalamus        ./ allmeta.braintotal
allmeta.corpus_callosum_normed = allmeta.corpus_callosum ./ allmeta.braintotal
allmeta.limbic_normed          = allmeta.limbic          ./ allmeta.braintotal
allmeta.subcortex_normed       = allmeta.subcortex       ./ allmeta.braintotal
allmeta.neocortex_normed       = allmeta.neocortex       ./ allmeta.braintotal
allmeta.cerebellum_normed      = allmeta.cerebellum      ./ allmeta.braintotal

### Metadata subgroups

ukidsmeta = view(allmeta, ukids, :)
oldkidsmeta = view(allmeta, oldkids, :)

@warn "Done!"
