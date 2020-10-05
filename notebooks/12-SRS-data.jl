using DrWatson; @quickactivate "ResonancePaper"

using ECHOAnalysis
using DataFrames
using CSV
using TOML: parsefile

srs =  CSV.File(datadir("metadata", "SRS_Data_All.csv")) |> DataFrame
rename!(srs, :studyID=>:subject)
# needs to have ENV["AIRTABLE_KEY"] = <key>
samplemeta = airtable_metadata()

# merge with subject metadata

allmeta = leftjoin(unique(samplemeta), srs, on=[:subject,:timepoint])

count(row-> !any(ismissing, row[["batch_16S", "SchoolageSRS::timepoint"]]), eachrow(allmeta))
count(row-> !any(ismissing, row[["batch_16S", "PreschoolSRS::timepoint"]]), eachrow(allmeta))

count(row-> !any(ismissing, row[["batch", "SchoolageSRS::timepoint"]]), eachrow(allmeta))
count(row-> !any(ismissing, row[["batch", "PreschoolSRS::timepoint"]]), eachrow(allmeta))
