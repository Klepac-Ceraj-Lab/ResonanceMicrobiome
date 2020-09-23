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

using DrWatson
@quickactivate "ResonancePaper"

using ECHOAnalysis
using TOML: parsefile
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

# needs to have ENV["AIRTABLE_KEY"] = <key>
samplemeta = airtable_metadata()

# merge with subject metadata

allmetadata = leftjoin(unique(samplemeta), subjectmeta, on=[:subject,:timepoint])
allmetadata.cogAssessment = [(ismissing(x) || x == "None") ? missing : x for x in allmetadata.cogAssessment]

# ## Brain Data
#
# We also have tables of brain volumes for many of our subjects.

brainfiles = config["tables"]["brain_structure"]

# Freesurfer is a way of doing segmentation
freesurfer = CSV.File(brainfiles["brainvolumes"]) |> DataFrame

# fix subjectID
function fixfreesurfersubject!(table)
    table.subject = map(table.ID) do id
        m = match(r"^sub-BAMBAM(\d+)$", id)
        isnothing(m) && error(id)
        parse(Int, m.captures[1])
    end
    select!(table, Not(:ID))
    return table
end

fixfreesurfersubject!(freesurfer)

for (n, col) in pairs(eachcol(freesurfer))
    eltype(col) <: Number && continue
    newcol = Union{Float64, Missing}[]
    for e in col
        if ismissing(e) || e == "#REF!"
            push!(newcol, missing)
        else
            push!(newcol, parse(Float64, e))
        end
    end
    freesurfer[!,n] = newcol
end

rename!(freesurfer, [   
    "Cerebral Spinal Fluid" => "csf",
    "Total Grey Matter Volume" => "gray_matter",
    "Cortical White Matter Volume" => "white_matter",
    "Brain-Stem" => "brainstem",
    "Total Intracranial Volume"=> "braintotal",
    "corpus callosum" => "corpus_callosum",
    ])

rename!(lowercase, freesurfer)

fs_keep = [
    "subject",
    "timepoint",
    "braintotal",
    "white_matter",
    "gray_matter",
    "csf",
    "brainstem",
    "hippocampus",
    "thalamus",
    "corpus_callosum",
    "limbic",
    "subcortex",
    "neocortex",
    "cerebellum"
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
    
