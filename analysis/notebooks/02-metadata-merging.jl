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

# Subject Metadata is stored in a CSV and can be easily loaded

using CSV
using DataFrames

include("airtable.key")
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

allmetadata = join(samplemeta, subjectmeta, on=[:subject,:timepoint], kind=:left)


# ## Brain Data
#
# We also have tables of brain volumes for many of our subjects.

brainfiles = config["tables"]["brain_structure"]
brainvol = CSV.read(brainfiles["lowres"])

# remove spaces from columns names
rename!(brainvol, map(names(brainvol)) do n
                        replace(String(n), " "=>"_") |> lowercase |> Symbol
                    end)
rename!(brainvol, :study_id => :subject)

# We need to fix the subjects - the letters represent timepoints -
# using the `resolve_letter_timepoint` function.

## convert letter timepoint into number
brainsid = resolve_letter_timepoint.(brainvol.subject)

brainvol.subject = subject.(brainsid)
brainvol.timepoint = timepoint.(brainsid)

allmetadata = join(allmetadata, brainvol, on=[:subject,:timepoint], kind=:left)


# And now the same thing for the high resolution scan table:

hires = CSV.read(brainfiles["hires"])
hires2 = CSV.read(brainfiles["hires2"])
rename!(hires, Dict(:ID=>:subject, :Timepoint=>:timepoint))
rename!(hires2, Dict(:ID=>:subject, :Timepoint=>:timepoint))


hr2_samples = resolve_letter_timepoint.(string.(hires2.subject))
hires2.subject = subject.(hr2_samples)
hires2.timepoint = timepoint.(hr2_samples)

# don't want to replicate :age column
select!(hires, Not(:Age))
@assert names(hires) == names(hires2)
hires = vcat(hires, hires2)
unique!(hires, [:subject,:timepoint])

# There are a lot of individual brain regions that are separated in this table,
# and the right and left hemispheres are distinguished.
# For the most part, we're not going to need this level of specificity,
# but we can group individual brain regions
# and combine left / right hemispheres.
# I'll also make a column with the total brain volume for later normalization.

mapping = CSV.read(brainfiles["hires_key"])
hires.hires_total = [sum(row[3:end]) for row in eachrow(hires)]

cols_seen = let ns = lowercase.(String.(names(hires)))
    cols_seen = Int[]
    by(mapping, :region) do region
        fs = lowercase.(region.feature)
        cols = findall(n-> any(f-> occursin(f, n), fs), ns)
        append!(cols_seen, cols)
        hires[!, Symbol(first(region.region))] = [sum(row[cols]) for row in eachrow(hires)]
        true # need something for the `by`
    end
    cols_seen
end

## join with other metadata
allmetadata = join(allmetadata, hires, on=[:subject,:timepoint], kind=:left)

## Write for easy referemce
CSV.write(config["tables"]["joined_metadata"], allmetadata)
