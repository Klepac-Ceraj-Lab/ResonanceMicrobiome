using DrWatson
@quickactivate "ResonancePaper"
using ECHOAnalysis
using DataFrames
using CSV
using TOML: parsefile
using AlgebraOfGraphics, AbstractPlotting

srs =  CSV.File(datadir("metadata", "SRS_Data_All.csv")) |> DataFrame
rename!(srs, :studyID=>:subject)

# change "" to actual key if environmental variable not set
ENV["AIRTABLE_KEY"] = get(ENV, "AIRTABLE_KEY", "")
samplemeta = airtable_metadata()
unique!(samplemeta, [:subject, :timepoint])
allmeta = CSV.File(datadir("metadata", "joined.csv")) |> DataFrame

# merge with subject metadata

allmeta = leftjoin(allmeta, srs, on=[:subject,:timepoint])
allmeta = leftjoin(allmeta, samplemeta, on=[:subject,:timepoint], makeunique=true)

names(allmeta)[findall(n-> occursin("16", n), names(allmeta))]

@info "16S has PreschoolSRS"
count(row-> !any(ismissing, row[["batch_16S", "PreschoolSRS::timepoint"]]), eachrow(allmeta)) |> println
@info "16S has SchoolageSRS"
count(row-> !any(ismissing, row[["batch_16S", "SchoolageSRS::timepoint"]]), eachrow(allmeta)) |> println

@info "mgx has PreschoolSRS"
count(row-> !any(ismissing, row[["batch", "PreschoolSRS::timepoint"]]), eachrow(allmeta)) |> println
@info "mgx has SchoolageSRS"
count(row-> !any(ismissing, row[["batch", "SchoolageSRS::timepoint"]]), eachrow(allmeta)) |> println

hassrs = filter(row-> !all(ismissing, row[["PreschoolSRS::timepoint", "SchoolageSRS::timepoint"]]), allmeta)
hassrs.correctedAgeYears = hassrs.correctedAgeDays ./ 365


hist = data(hassrs) * mapping(:correctedAgeYears) * AlgebraOfGraphics.histogram |> draw



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
