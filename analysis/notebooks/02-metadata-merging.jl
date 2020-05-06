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
config = parsefile("data/data.toml")
for (key, value) in config
    println(key,":")
    println("\t",value)
end

##
# Subject Metadata is stored in a CSV and can be easily loaded

using CSV
using DataFrames
using PrettyTables

# make a formatter for PrettyTables to round long numbers
rounder = Dict(0 => (v,i) -> typeof(v) <: AbstractFloat ? round(v,digits=3) : v)
# make a row filter that picks ~15 rows
randrowfilter(data, i) = rand() < (1 / size(data, 1)) * 15

@ptconfclean # clear previous configuration
@ptconf formatter = rounder nosubheader=true screen_size=(20,120) filters_row=(randrowfilter,)
# set configuration for printing tables
@ptconf formatter = rounder nosubheader=true screen_size=(20,120)

##

include("airtable.key")
subjectmeta = echo_
# pretty print table
@pt allmeta

## Sample metadata

In addition to the FilemakerPro database,
we also have metdata info stored for each of the samples that are processed.
In this case, `timepoint` and `subject` do not uniquely identify samples,
since we can have multiple samples per timepoint.
`sample` IDs should be unique though.

```julia; results="hidden"
samples = CSV.read(config["tables"]["fecal_samples"]["path"])
rename!(samples, [:TimePoint=>:timepoint, :DOC=>:date, :SubjectID=>:subject, :SampleID=>:sample])

# convert to longform
samples = melt(samples, [:subject, :timepoint, :sample], variable_name=:metadatum)
dropmissing!(samples, [:subject, :sample])
samples[!, :parent_table] .= "FecalProcessing"
```

Finally, we want to make sure the maternal samples
also have a value for `ageLabel`.

```julia
momlabels = by(samples, [:subject, :timepoint]) do sample
    if any(s-> startswith(s, "M"), skipmissing(sample.sample))
        outdf = filter(row-> !ismissing(row.sample) && startswith(row.sample, "M"), sample)
        outdf.metadatum .= :ageLabel
        outdf.value .= "mom"
        outdf.parent_table .= "Calculated"
    else
        outdf = DataFrame(first(sample))

    end
    return unique(outdf[:, [:metadatum, :value, :parent_table, :sample]])
end
filter!(row-> row.value == "mom", momlabels)
@pt momlabels[.!ismissing.(momlabels.value), :]
```
```julia; results="hidden"
# add to main sample dataframe (while reordering columns)
append!(samples, momlabels[!, names(samples)])
```

## Brain Data

We also have tables of brain volumes for many of our subjects.

```julia
brainfiles = config["tables"]["brain_structure"]
brainvol = CSV.read(brainfiles["lowres"])

# remove spaces from columns names
names!(brainvol, map(names(brainvol)) do n
                        replace(String(n), " "=>"_") |> lowercase |> Symbol
                    end)
rename!(brainvol, :study_id => :subject)

# Convert to longform
brainvol = stack(brainvol, [:white_matter_volume, :grey_matter_volume, :csf_volume], :subject, variable_name=:metadatum)
@pt brainvol
```

We need to fix the subjects - the letters represent timepoints -
using the `resolve_letter_timepoint` function.

```julia
# convert letter timepoint into number
brainsid = resolve_letter_timepoint.(brainvol.subject)
brainsid[1:5]
```

```julia
brainvol.subject = subject.(brainsid)
brainvol.timepoint = timepoint.(brainsid)
brainvol.sample = sampleid.(brainsid)
# This uses slightly different syntax, because the column doesn't exist.
# The `.=` broadcast assigment makes every row in the column the same string
brainvol[!, :parent_table] .= "brainVolume"

@pt brainvol
```

And now the same thing for the high resolution scan table:

```julia; results="hidden"
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

@pt hires
```

There are a lot of individual brain regions that are separated in this table,
and the right and left hemispheres are distinguished.
For the most part, we're not going to need this level of specificity,
but we can group individual brain regions
and combine left / right hemispheres.
I'll also make a column with the total brain volume for later normalization.

```julia
mapping = CSV.read(brainfiles["hires_key"])
@pt mapping
```
```julia
hires.hires_total = [sum(row[3:end]) for row in eachrow(hires)]

cols_seen = let ns = lowercase.(String.(names(hires)))
    cols_seen = Int[]
    by(mapping, :region) do region
        fs = lowercase.(region.feature)
        cols = findall(n-> any(f-> occursin(f, n), fs), ns)
        append!(cols_seen, cols)
        hires[!, Symbol(first(region.region))] = [sum(row[cols]) for row in eachrow(hires)]
    end
    cols_seen
end
```



```julia
# convert to longform
hires = melt(hires, [:subject, :timepoint], variable_name=:metadatum)
hires[!, :parent_table] .= "hiresVolumes"
```

We can only concatenate tables if they all have the same columns,
so I'll add a `sample` ID to all of the other observations
to match what's in the `samples` DataFrame.
The fecal sample `sample` IDs are build from the `subject` ID and `timepoint`,
so I'll do the same for other observations.

_Note_: Fecal sample `sample` IDs are built as follows...

```
C0596_1F_1A
C = Child (or M = Mother)
0596 = SubjectID
1 = TimePoint (converted from A=1, B=2, etc.)
F = Fecal i.e. Genotek sample (or E = ethanol sample)
1 = CollectionRep (if multiple samples for same TimePoint)
A = AliquotRep (each fecal Genotek sample is aliquoted into 2 or 4 smaller cryovials)
    A,B or A,B,C,D ; SOP is to process AliquotRep A for DNA extractions
```

```julia; results="hidden"
allmeta[!,:sample] = map(r->
        "C" * lpad(string(r[:subject]), 4, "0") * "_$(Int(floor(r[:timepoint])))M",
        eachrow(allmeta))

brainvol[!,:sample] = map(r->
        "C" * lpad(string(r[:subject]), 4, "0") * "_$(Int(floor(r[:timepoint])))M",
        eachrow(brainvol))

hires[!,:sample] = map(r->
        "C" * lpad(string(r[:subject]), 4, "0") * "_$(Int(floor(r[:timepoint])))M",
        eachrow(hires))

```

And then concatenate all the tables together

```julia
allmeta = vcat(allmeta, brainvol, hires, samples)
# reorder columns
permutecols!(allmeta, [:sample, :subject, :timepoint, :metadatum, :value, :parent_table])
# remove rows with missing values
dropmissing!(allmeta)
allmeta.metadatum = string.(allmeta.metadatum)
allmeta.value = string.(allmeta.value)
@pt allmeta
```

## Save to SQLite

```julia
SQLite.drop!(metadb, "allmetadata", ifexists=true)
SQLite.dropindex!(metadb, "allmetadata_subject_idx", ifexists=true)
SQLite.dropindex!(metadb, "allmetadata_metadatum_idx", ifexists=true)

allmeta |> SQLite.load!(metadb, "allmetadata")
SQLite.createindex!(metadb, "allmetadata", "allmetadata_subject_idx", "subject", unique=false)
SQLite.createindex!(metadb, "allmetadata", "allmetadata_metadatum_idx", "metadatum", unique=false)
```

## Converting to wide-form

The `ECHOAnalysis` repo has a function
to convert long-form data to wide-form,
where each sample is on a singe row.
For metadata linked to a timpoint,
values are only applied to the sample linked to that timepoint.
For metadata not linked to a timepoint,
the value is applied to every timepoint for that subject.


```julia
widemeta = widemetadata(allmeta, unique(samples.sample))
@pt widemeta
```

Unfortunately, `sqlite` doesn't handle the really wide table well.
But it's not very big,
So I'll just store it as a regular CSV file.

```julia
CSV.write(config["tables"]["sample_metadata_wide"]["path"], widemeta)
```
