# ---
# title: "Notebook 4: Quality control of metagenomes"
# author: "Kevin Bonham, PhD"
# options:
    # line_width : 120
    # wrap : false
# ---
# 
# All of the metagenomes were processed
# using tools from the [bioBakery](https://bitbucket.org/biobakery/biobakery/wiki/Home).
# 
# Note - these files needed for this analysis
# are not included in the zenodo data repository,
# but are avaliable on request.

ENV["GKSwstype"] = "100"
using ECHOAnalysis
using DataFrames
using SQLite
using StatsMakie
using AbstractPlotting.MakieLayout
using PrettyTables
using CSV
using Pkg.TOML: parsefile

rounder = Dict(0 => (v,i) -> typeof(v) <: AbstractFloat ? round(v,digits=3) : v)
# print ~15 random rows
randrowfilter(data, i) = rand() < (1 / size(data, 1)) * 15
@ptconfclean # clear previous configuration
@ptconf formatter = rounder nosubheader=true screen_size=(20,120) filters_row=(randrowfilter,)

config = parsefile("data/data.toml")
widemeta = ECHOAnalysis.getmgxmetadata()
figures = config["output"]["figures"]
tables = config["output"]["tables"]
isdir(figures) || mkpath(figures)
isdir(tables) || mkpath(tables)

# ## Quality control

# First, I'll look at the QC results from `kneaddata`.

qc_files = let qcfiles=[]
    for (root, dirs, files) in walkdir(config["files"]["biobakery"]["path"])
        occursin("kneaddata", root) || continue
        filter!(files) do f
            occursin("read_counts", f)
        end
        append!(qcfiles, joinpath.(root, files))
    end
    qcfiles
end

## make a DataFrame of the qc table from the first batch
qc = CSV.read(qc_files[1])
## sum paired end reads
qc[!,:raw] = qc[!,Symbol("raw pair1")] .+ qc[!,Symbol("raw pair2")]
qc[!,:trimmed] = qc[!,Symbol("trimmed pair1")] .+ qc[!,Symbol("trimmed pair2")]
qc[!,:final] = qc[!,Symbol("final pair1")] .+ qc[!,Symbol("final pair2")]
qc[!, :batch] .= "batch001"
qc = qc[!, [:Sample, :raw, :trimmed, :final, :batch]]

## loop through the rest and concatenate to the first one
for f in qc_files[2:end]
    df = CSV.read(f)
    ## get batch name from file
    df[!, :batch] .= match(r"(batch\d+)", f).captures[1]
    ## column names are different depending on whether reads were concatenated or paired
    if any(x-> occursin("pair", x), String.(names(df)))
        df[!,:raw] = df[!,Symbol("raw pair1")] .+ df[!,Symbol("raw pair2")]
        df[!,:trimmed] = df[!,Symbol("trimmed pair1")] .+ df[!,Symbol("trimmed pair2")]
        df[!,:final] = df[!,Symbol("final pair1")] .+ df[!,Symbol("final pair2")]
        df = df[!, [:Sample, :raw, :trimmed, :final, :batch]]
    else
        rename!(df, [Symbol("raw single")=>:raw,
                    Symbol("trimmed single")=>:trimmed,
                    Symbol("final single")=>:final])
        df = df[!, [:Sample, :raw, :trimmed, :final, :batch]]
    end
    global qc = vcat(qc, df)
end

@pt qc

# To keep the formatting of sample IDs consistant across data types,
# I'll use the `stoolsample` function.

let samples = qc[!, :Sample]
    for i in eachindex(samples)
        s = replace(samples[i], "_kneaddata"=> "")
        samples[i] = sampleid(stoolsample(s))
    end
end

@pt qc

# Now let's take a look at them with some plots.

## sort by batch, then by raw read count
sort!(qc, [:batch, :raw])

# TODO: Switch to Makie
bar(x=qc.Sample, hcat(qc.raw, qc.final),
    xaxis="Samples", yaxis= "Count", legend=:topleft,
    title = "QC from Kneaddata", label=["Raw" "Final"],
    linecolor=:match)


savefig(joinpath(figures, "04-knead-qc.png"))

# These are a little more variable than I'd like.
# Let's take a look at their properties:

using Statistics

qc_stats = by(qc, :batch) do df
                DataFrame(
                  mean = round(mean(df.final) / 1e6, digits=2),
                  med  = round(median(df.final) / 1e6, digits=2),
                  max  = round(maximum(df.final) / 1e6, digits=2),
                  min  = round(minimum(df.final) / 1e6, digits=2),
                  )
end
CSV.write(joinpath(tables, "qc_stats.csv"), qc_stats)
@pt qc_stats

# According to Andre Comeau of [Integrated Microbiome Research](http://www.imr.bio)
# (where our sequencing is done):

# > However, even with sample normalization, the best value that kits/protocols can obtain is about a 2-fold difference from the mean...which then means if your average # of reads is, for example, 8 M, then you'll see samples up to 2-fold higher and 2-fold lower at max = about 4-16 M range.
# >
# > Now added on top of that is that each NextSeq run is independent and the loading (cluster density, which is a bit of an art) tends to vary a little bit, so overall output per sample also varies there too, but usually within a tighter "true" 2-fold range. Hence this explains a bunch of the variation.
