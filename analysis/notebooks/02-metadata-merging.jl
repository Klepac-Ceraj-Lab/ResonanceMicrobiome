# ---
# title: "Notebook 2: Working with Metadata"
# author: "Kevin Bonham, PhD"
# options:
#     line_width : 120
#     wrap : false
# ---
#

# ## Accessing TOML data in julia
#
# Information about the locations of data are found in `Data.toml`.
# If you've downloaded the data from Zenodo,
# be sure to update the paths in that file for this code to run successfully.
# Parsing this file gives a set of nested key:value pairs.
#
# Extra code for much of this analysis is found in the `ECHOAnalysis` julia package.
# The docs can be [found here](https://klepac-ceraj-lab.github.io/echo_analysis/dev/).

using ECHOAnalysis
using Pkg.TOML: parsefile
config = parsefile("Data.toml")

for (key, value) in config
    println(key,":")
    println("\t",value)
end

## Subject Metadata is stored in a CSV and can be easily loaded

using CSV
using DataFrames

subjectmeta = CSV.File(config["tables"]["subject_metadata"]) |> DataFrame
rename!(subjectmeta, :subjectID=> :subject)

# ## Sample metadata
#
# In addition to the FilemakerPro database,
# we also have metdata info stored for each of the samples that are processed.
# These can be loaded directly from the airtable database
# if you have the API key.
# If you dowloaded this table from Zenodo, skip this step

include("airtable.key"); # ENV["AIRTABLE_KEY"] = <key>
samplemeta = airtable_metadata()

# merge with subject metadata

allmetadata = leftjoin(unique(samplemeta), subjectmeta, on=[:subject,:timepoint])
allmetadata.cogAssessment = [(ismissing(x) || x == "None") ? missing : x for x in allmetadata.cogAssessment]

# ## Brain Data
#
# We also have tables of brain volumes for many of our subjects.

brainfiles = config["tables"]["brain_structure"]

# Freesurfer is another way of doing segmentation
# We have 2 different versions of freesurfer tables to merge


freesurfer5 = CSV.File(brainfiles["freesurferv5"]) |> DataFrame
freesurfer6 = CSV.File(brainfiles["freesurferv6"]) |> DataFrame

# Some sanity checks
@warn "Checking freesurfer tables"

@info "v5 - v6 cols" setdiff(names(freesurfer5), names(freesurfer6))
@info "v6 - v5 cols" setdiff(names(freesurfer6), names(freesurfer5))

@assert names(freesurfer5) == names(freesurfer6)

# fix subjectID
function fixfreesurfersubject!(table)
    table.subject = map(table.ID) do id
        m = match(r"^sub-BAMBAM(\d+)$", id)
        isnothing(m) && error(id)
        parse(Int, m.captures[1])
    end
    return table
end

fixfreesurfersubject!(freesurfer5)
fixfreesurfersubject!(freesurfer6)

let ns = names(freesurfer5)
    idx_all = findall(n-> eltype(freesurfer5[!,n]) <: Real && all(==(0.), freesurfer5[!,n]), ns)
    idx_any = findall(n-> eltype(freesurfer5[!,n]) <: Real && any(==(0.), freesurfer5[!,n]), ns)
    @info "v5 all zeros" ns[idx_all]
    @info "v5 any zeros" setdiff(ns[idx_any], ns[idx_all])
end

let ns = names(freesurfer6)
    idx_all = findall(n-> eltype(freesurfer6[!,n]) <: Real && all(==(0.), freesurfer6[!,n]), ns)
    idx_any = findall(n-> eltype(freesurfer6[!,n]) <: Real && any(==(0.), freesurfer6[!,n]), ns)
    @info "v6 all zeros" ns[idx_all]
    @info "v6 any zeros" setdiff(ns[idx_any], ns[idx_all])
end

let dupes = intersect(collect(zip(freesurfer5.subject, freesurfer5.timepoint)),collect(zip(freesurfer6.subject, freesurfer6.timepoint)))
    fs5 = filter([:subject, :timepoint] => (s,t) -> (s,t) ∈ dupes, freesurfer5)
    fs6 = filter([:subject, :timepoint] => (s,t) -> (s,t) ∈ dupes, freesurfer6)

    ns = names(fs5)
    @assert ns == names(fs6)
    @assert size(fs5) == size(fs6)
    for (r5, r6) in zip(eachrow(fs5), eachrow(fs6))
        @info "version mismatches"   ns[findall(n-> r5[n] != r6[n], ns)]
    end
    fs5, fs6
end

freesurfer = vcat(freesurfer5, freesurfer6)
unique!(freesurfer, [:subject, :timepoint])

freesurfer.hippocampus = freesurfer."Left-Hippocampus" .+ freesurfer."Right-Hippocampus"
freesurfer.caudate = freesurfer."Left-Caudate" .+ freesurfer."Right-Caudate"
freesurfer.putamen = freesurfer."Left-Putamen" .+ freesurfer."Right-Putamen"
freesurfer.pallidum = freesurfer."Left-Pallidum" .+ freesurfer."Right-Pallidum"
freesurfer.thalamus = freesurfer."Left-Thalamus-Proper" .+ freesurfer."Right-Thalamus-Proper"
freesurfer.amygdala = freesurfer."Left-Amygdala" .+ freesurfer."Right-Amygdala"
freesurfer.corpus_callosum = freesurfer.CC_Posterior .+ freesurfer.CC_Mid_Posterior .+ freesurfer.CC_Central .+ freesurfer.CC_Mid_Anterior .+ freesurfer.CC_Anterior

rename!(freesurfer, [
    "CSF" => "csf",
    "TotalGrayVol" => "gray_matter",
    "Brain-Stem" => "brainstem",
    "BrainSegVol"=> "braintotal",
    "CorticalWhiteMatterVol" => "white_matter"
    ])

fs_keep = [
    "subject",
    "timepoint",
    "braintotal",
    "white_matter",
    "gray_matter",
    "csf",
    "brainstem",
    "hippocampus",
    "caudate",
    "putamen",
    "pallidum",
    "thalamus",
    "amygdala",
    "corpus_callosum"
]

select!(freesurfer, fs_keep)
allmetadata = join(allmetadata, freesurfer, on=[:subject,:timepoint], kind=:left)


## Write for easy referemce
CSV.write(config["tables"]["joined_metadata"], allmetadata)
CSV.write("/home/kevin/Desktop/hasstool.csv", unique(samplemeta[map(s-> startswith(s, "C"), samplemeta.sample), [:subject, :timepoint]]))

ukids, oldkids = let samples = Set(sampleid.(uniquetimepoints(allmetadata.sample, takefirst=true, samplefilter=iskid)))
    (map(row-> !ismissing(row.ageLabel) && in(row.sample, samples), eachrow(allmetadata)),
    map(row-> !ismissing(row.ageLabel) && in(row.sample, samples) && row.ageLabel != "1 and under", eachrow(allmetadata)))
end

noreps = let samples = Set(sampleid.(uniquetimepoints(stoolsample.(allmetadata.sample), takefirst=false, samplefilter=iskid)))
    map(row-> !ismissing(row.ageLabel) && in(row.sample, samples), eachrow(allmetadata))
end

norepsmeta = view(allmetadata, noreps, :)
ukidsmeta = view(allmetadata, ukids, :)
oldkidsmeta = view(allmetadata, oldkids, :)

print("All: ","\n\t",
    "N samples: ", size(norepsmeta, 1), "\n\t",
    "Has scan:  ", count(row-> !ismissing(row.braintotal), eachrow(norepsmeta)), "\n\t",
    "Has cog:   ", count(row-> !ismissing(row.cogScore), eachrow(norepsmeta)), "\n",
    "Unique: ","\n\t",
    "N samples: ", size(ukidsmeta, 1), "\n\t",
    "Has scan:  ", count(row-> !ismissing(row.braintotal), eachrow(ukidsmeta)), "\n\t",
    "Has cog:   ", count(row-> !ismissing(row.cogScore), eachrow(ukidsmeta)), "\n",
    "> 1 yo: ","\n\t",
    "N samples: ", size(oldkidsmeta, 1), "\n\t",
    "Has scan:  ", count(row-> !ismissing(row.braintotal), eachrow(oldkidsmeta)), "\n\t",
    "Has cog:   ", count(row-> !ismissing(row.cogScore), eachrow(oldkidsmeta)), "\n"
    )
    
