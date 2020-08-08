# ---
# title: "Notebook 3: Metadata Analysis"
# author: "Kevin Bonham, PhD"
# options:
#     line_width : 120
#     wrap : false
# ---

# ## Getting Data

using ECHOAnalysis
using DataFrames
using StatsMakie
using AbstractPlotting
using AbstractPlotting.MakieLayout
using PrettyTables
using CSV
using Pkg.TOML: parsefile
rounder = Dict(0 => (v,i) -> typeof(v) <: AbstractFloat ? round(v,digits=3) : v)

randrowfilter(data, i) = rand() < (1 / size(data, 1)) * 15
@ptconfclean # clear previous configuration
@ptconf formatter = rounder nosubheader=true screen_size=(20,120) filters_row=(randrowfilter,)

config = parsefile("data/data.toml")
widemeta = ECHOAnalysis.getmgxmetadata()
figures = config["output"]["figures"]
isdir(figures) || mkpath(figures)

# As an example, how many unique subjects
# do we have any metadata for?

# the |> is pipe syntax, the following is the same as
# `length(unique(widemeta[:subject]))`

unique(widemeta.subject) |> length

# (**Note**: We have metadata for more subjects,
# but this table was created
# for only the subjects that have provided fecal samples)

# How many samples for each subject?

sampleinfo = by(widemeta, :subject) do df
    (;nsamples = size(df,1))
end

# TODO: change to Makie
histogram(sampleinfo.nsamples, legend=false,
    title="Samples per Subject ID",
    xlabel="# of fecal samples", ylabel="# of subjects",
    xticks=1:12)

savefig(joinpath(figures, "03-samples-per-subject.svg"))

# Wow - there are a couple of subjects that have a lot of samples.
# Which subjects are those?

highsamplers = filter(row-> row.nsamples > 4, sampleinfo).subject

# To see what samples are from those folks that gave a bunch:

highsamplersinfo = filter(widemeta) do row
    # find rows from subjects in the highly sampled pool
    in(row.subject, highsamplers)
end

sort!(highsamplersinfo, [:subject, :timepoint])

@pt highsamplersinfo[!, [:sample, :subject, :timepoint]]

# So a bunch of these are where
# multiple samples were given for the same timepoint (eg `C0016_3F_1A` and `_2A`)
# and/or both genotek (`F`) and enthanol (`E`) samples.

# _Note_: after `batch006`, SOP is to send 1 `ALiquotRep` of 1 `CollectionRep` for
# each `timepoint` (i.e. only send `C0202_4F_1A` for mgx sequencing, not `C0202_4F_2A`
# or `C0202_4F_1B` or `C0202_4E_1A`) unless otherwise noted.

# ### Unique samples

# Except for later quality control,
# we don't actually want to analyze replicates.
# And we want to focus on only samples collected in genotek tubes,
# not ethanol (so ones with an F in the second ID slot).

# Using the `stoolsample` function,
# I can get just the relevant info for filtering on the first sample.

samples = stoolsample.(widemeta.sample)
samples[1:4]
seen = Timepoint[]
usamples = StoolSample[]

sort!(samples)

map(samples) do s
    # skip ethanol samples
    sampletype(s) == "ethanol" && return nothing
    subtp = Timepoint("id", subject(s), timepoint(s))
    if !in(subtp, seen)
        push!(seen, subtp)
        push!(usamples, s)
    end
end

# The `ECHOAnalysis` module also has [a function](https://klepac-ceraj-lab.github.io/echo_analysis/dev/metadata_handling/#ECHOAnalysis.uniquesamples-Tuple{AbstractArray{#s17,1}%20where%20#s17%3C:NamedTuple})
# that does this, and has a bunch of other options too.

usamples2 = uniquetimepoints(samples)
## this would throw an error if false
@assert usamples2 == usamples

# To get just the metadata for these uniqe samples:

umeta = let us = Set(sampleid.(usamples))
    filter(row-> row.sample in us, widemeta)
end
sort!(umeta, :sample)


# Metadata in this wide form is useful for things like plotting.
# For example, let's look at the age-adjusted cognitive scores.
# First, we'll filter on kids that have cognitive scores and ages:

using ColorSchemes
p1 = ColorSchemes.Set1_8.colors

toplot = filter(row-> !any(ismissing, [row.cogScore, row.correctedAgeDays]), umeta)
scatter(toplot.correctedAgeDays ./ 365, toplot.cogScore, group=toplot.cogAssessment,
    color=p1[2:end]', xlabel="age (years)", ylabel="Composite Score",
    title="Cognitive Assessments", legend=:topright,
    )
savefig(joinpath(figures, "03-cogscore_age_scatter.svg"))

using Statistics
m = mean(toplot.cogScore)
s = std(toplot.cogScore)

toplot.zscore = map(x-> (x - m)/s, toplot.cogScore)

scatter(toplot.correctedAgeDays ./ 365, toplot.zscore,
    group = toplot.cogAssessment, color=p1[2:end]',
    ylabel = "zscore", xlabel = "age", legend=:topright,
    title="Cognitive Assessments")
savefig(joinpath(figures, "03-cogscore_zscore_age_scatter.svg"))

Checking if any of the outliers are correlated with batch number:

toplot = filter(row-> !any(ismissing, [row.Mgx_batch, row.correctedAgeDays, row.cogScore]), umeta)
toplot.Mgx_batch
scatter(toplot.correctedAgeDays ./ 365, toplot.cogScore, primary=false,
    zcolor=toplot.Mgx_batch, xlabel="age (years)", ylabel="Composite Score", group=toplot.Mgx_batch,
    title="Cognitive Assessments", legend=:topright)
savefig(joinpath(figures, "03-cogscore_age_scatter_batchnum.svg"))
