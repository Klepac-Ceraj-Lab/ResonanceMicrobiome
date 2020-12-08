# Want to get samples for which we have mother / child pairs
# Targeting ~ 250 samples

using DataFrames
using Airtable
using Chain
using CSV

# replace "" with actual key if it's not in ENV
get(ENV, "AIRTABLE_KEY", "")

batch2string(b) = ismissing(b) ? missing : parse(Int, match(r"Batch (\d+)", b).captures[1])
# this will eventually come from ECHOAnalysis.jl when I work out dep problems
function airtable_metadata(key=Airtable.Credential())
    df = DataFrame()
    req = Airtable.get(key, "/v0/appyRaPsZ5RsY4A1h", "Master"; view="ALL_NO_EDIT")
    foreach(req.records) do record
        f = record.fields
        push!(df, f, cols=:union)
    end
    while haskey(req, :offset)
        @info "Making another request"
        req = Airtable.get(key, "/v0/appyRaPsZ5RsY4A1h/", "Master"; view="ALL_NO_EDIT", offset=req.offset)
        foreach(req.records) do record
            f = record.fields
            push!(df, f, cols=:union)
        end
        sleep(0.250)
    end

    rename!(df, :SampleID=>:sample)
    df.subject = tryparse.(Int, df.SubjectID)
    df.timepoint = tryparse.(Int, df.TimePoint)

    df.mgxbatch = map(batch2string, df.Mgx_batch)
    df."16Sbatch" = map(batch2string, df."16S_batch")

    return df[!, Not(["SubjectID", "TimePoint", "Mgx_batch", "16S_batch"])]
end

df = airtable_metadata()
filter!(:Fecal_EtOH=> ==("E"), df)

@chain df begin
    groupby(:subject)
    transform!(nrow => :subjectcount)
    filter!(:subjectcount => >(1), _)
    sort!([:subjectcount, :subject])
end

@chain df begin
    groupby(:subject)
    transform!(:Mother_Child => (mc-> length(unique(mc)) == 2) => :hasboth)
    filter!(:hasboth => Bool, _)
    filter!(:CollectionRep => ==("1"), _)
end


CSV.write("/home/kevin/Desktop/metabolomics_samples.csv", select!(df, [:sample, :subject, :timepoint, :StorageBox, :EthanolStorage_old, :subjectcount]))
