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

allmetadata = join(unique(samplemeta), subjectmeta, on=[:subject,:timepoint], kind=:left)


# ## Brain Data
#
# We also have tables of brain volumes for many of our subjects.

brainfiles = config["tables"]["brain_structure"]

# Freesurfer is another way of segmentation

freesurfer = CSV.File(brainfiles["freesurfer"]) |> DataFrame
# fix subjectID
freesurfer.subject = map(freesurfer.ID) do id
    m = match(r"^sub-BAMBAM(\d+)$", id)
    isnothing(m) && error(id)
    parse(Int, m.captures[1])
end

# one subject was not segmented properly for some reason
filter!(row-> row.subject != 767, freesurfer)

freesurfer.hippocampus = freesurfer."Left-Hippocampus" .+ freesurfer."Right-Hippocampus"
freesurfer.caudate = freesurfer."Left-Caudate" .+ freesurfer."Right-Caudate"
freesurfer.putamen = freesurfer."Left-Putamen" .+ freesurfer."Right-Putamen"
freesurfer.pallidum = freesurfer."Left-Pallidum" .+ freesurfer."Right-Pallidum"
freesurfer.thalamus = freesurfer."Left-Thalamus-Proper" .+ freesurfer."Right-Thalamus-Proper"
freesurfer.amygdala = freesurfer."Left-Amygdala" .+ freesurfer."Right-Amygdala"
freesurfer.corpus_callosum = freesurfer.CC_Posterior .+ freesurfer.CC_Mid_Posterior .+ freesurfer.CC_Central .+ freesurfer.CC_Mid_Anterior .+ freesurfer.CC_Anterior
# these are the same thing in different versions of freesurfer - should be 0 in one column
@assert all(row-> xor(row.CerebralWhiteMatterVol == 0, row.CorticalWhiteMatterVol == 0), eachrow(freesurfer))
freesurfer.white_matter = freesurfer.CerebralWhiteMatterVol .+ freesurfer.CorticalWhiteMatterVol


rename!(freesurfer, [
    "CSF" => "csf",
    "TotalGrayVol" => "gray_matter",
    "Brain-Stem" => "brainstem",
    "BrainSegVol"=> "braintotal"
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

ukidsmeta = view(allmetadata, ukids, :)
oldkidsmeta = view(allmetadata, oldkids, :)

println(size(allmetadata, 1))
println(size(ukidsmeta, 1))
println(size(oldkidsmeta, 1))
println(count(row-> !ismissing(row.braintotal), eachrow(allmetadata)))
println(count(row-> !ismissing(row.braintotal), eachrow(ukidsmeta)))
println(count(row-> !ismissing(row.braintotal), eachrow(oldkidsmeta)))
println(count(row-> !ismissing(row.cogScore), eachrow(allmetadata)))
println(count(row-> !ismissing(row.cogScore), eachrow(ukidsmeta)))
println(count(row-> !ismissing(row.cogScore), eachrow(oldkidsmeta)))

