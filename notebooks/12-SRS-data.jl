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

@info "16S has PreschoolSRS"
count(row-> !any(ismissing, row[["batch_16S", "PreschoolSRS::timepoint"]]), eachrow(allmeta)) |> println
@info "16S has SchoolageSRS"
count(row-> !any(ismissing, row[["batch_16S", "SchoolageSRS::timepoint"]]), eachrow(allmeta)) |> println

@info "mgx has PreschoolSRS"
count(row-> !any(ismissing, row[["batch", "PreschoolSRS::timepoint"]]), eachrow(allmeta)) |> println
@info "mgx has SchoolageSRS"
count(row-> !any(ismissing, row[["batch", "SchoolageSRS::timepoint"]]), eachrow(allmeta)) |> println


## CBCL

cbcl = CSV.File(datadir("metadata", "has_CBCL.csv")) |> DataFrame
allmeta = leftjoin(unique(samplemeta), cbcl, on=[:subject,:timepoint])

@info "16S has younger CBCL"
count(row-> !any(ismissing, row[["batch_16S", "cbclAge"]]) && row["cbclAge"] == "younger", eachrow(allmeta)) |> println
@info "16S has older CBCL"
count(row-> !any(ismissing, row[["batch_16S", "cbclAge"]]) && row["cbclAge"] == "older", eachrow(allmeta)) |> println

@info "mgx has younger CBCL"
count(row-> !any(ismissing, row[["batch", "cbclAge"]]) && row["cbclAge"] == "younger", eachrow(allmeta)) |> println
@info "mgx has older CBCL"
count(row-> !any(ismissing, row[["batch", "cbclAge"]]) && row["cbclAge"] == "older", eachrow(allmeta)) |> println
