# Want to get samples for which we have mother / child pairs
# Targeting ~ 250 samples

using DataFrames
using Airtable
using Chain
using CSV

# replace "" with actual key if it's not in ENV
key = Airtable.Credential(get(ENV, "AIRTABLE_KEY", ""))

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

df = airtable_metadata(key)
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

samplessent = Set(["M0652_1E_1A","M0647_1E_1A","M0647_2E_1A","M0676_1E_1A","M0682_1E_1A","M0711_1E_1A","M0713_1E_1A","M0723_1E_1A","M0754_1E_1A","M0652_3E_1A","M0753_1E_1A","M0774_1E_1A","M0742_2E_1A","M0713_2E_1A","M0792_1E_1A","M0711_2E_1A","M0752_2E_1A","M0774_2E_1A","M0794_1E_1A","M0753_2E_1A","M0834_1E_1A","M0713_3E_1A","M0828_1E_1A","M0769_2E_1A","M0839_1E_1A","M0863_1E_1A","M0752_3E_1A","M0851_1E_1A","M0769_3E_1A","C0723_3E_1A","M0794_2E_1A","M0886_1E_1A","M0880_1E_1A","M0893_1E_1A","M0753_3E_1A","C0834_4E_1A","M1059_1E_1A","C1415_2E_1A","C0880_6E_1A","C1344_4E_1A","C1088_3E_1A","C1218_5E_1A","M1415_3E_1A","C1398_3E_1A","C1322_4E_1A","C1246_5E_1A","C1401_3E_1A","C0774_5E_1A","C0864_8E_1A","M0908_1E_1A","C0754_2E_1A","M0863_2E_1A","M0880_2E_1A","M0839_3E_1A","M0915_1E_1A","C0652_4E_1A","M0888_2E_1A","M0855_1E_1A","M0932_1E_1A","M0904_2E_1A","M0893_2E_1A","M0886_2E_1A","M0908_2E_1A","M0863_3E_1A","C0682_3E_1A","C0711_4E_1A","M0932_2E_1A","C0742_3E_1A","C0647_4E_1A","M0891_3E_1A","M0996_1E_1A","C0713_4E_1A","M1007_1E_1A","M0886_3E_1A","C0752_4E_1A","M0980_2E_1A","M1015_1E_1A","C0754_3E_1A","M0958_2E_1A","M0961_2E_1A","M1115_1E_1A","C0514_2E_1A","C0915_3E_1A","M1050_1E_1A","C0880_4E_1A","M1087_2E_1A","M1160_1E_1A","M1060_2E_1A","C0855_4E_1A","C0863_4E_1A","M1099_2E_1A","M1061_3E_1A","M1167_1E_1A","C0769_5E_1A","M1162_1E_1A","M1068_2E_1A","M1088_2E_1A","M1162_2E_1A","M1159_2E_1A","C0754_4E_1A","M1177_2E_1A","M1159_3E_1A","M1109_3E_1A","C1059_3E_1A","M1207_2E_1A","M1220_2E_1A","C1160_3E_1A","M1220_1E_1A","C0891_4E_1A","M1195_2E_1A","M1115_3E_1A","C0723_4E_1A","M1227_1E_1A","M1162_3E_1A","M1234_1E_1A","M1155_3E_1A","M1246_1E_1A","C0932_4E_1A","C0904_5E_1A","M1195_3E_1A","M1177_3E_1A","M1234_2E_1A","M1271_1E_1A","M1234_3E_1A","C1060_3E_1A","C0958_5E_1A","C0996_4E_1A","C1050_4E_1A","C1109_4E_1A","M1314_2E_1A","M1326_2E_1A","M1322_2E_1A","M1060_1E_1A","M1015_2E_1A","M1087_1E_1A","C0839_4E_1A","M1008_3E_1A","C0888_4E_1A","C0958_4E_1A","M1159_1E_1A","M1050_3E_1A","M1135_2E_1A","M1109_2E_1A","M1197_1E_1A","M1115_2E_1A","M1155_2E_1A","C0908_4E_1A","C1037_4E_1A","C0839_5E_1A","C0864_5E_1A","C1007_3E_1A","C0961_4E_1A","C1068_3E_1A","C0914_5E_1A","C0753_4E_1A","M1218_2E_1A","M1219_2E_1A","M1008_2E_1A","M0996_2E_1A","M0961_3E_1A","M0958_3E_1A","C0711_5E_1A","M1068_1E_1A","M1061_2E_1A","C0851_4E_1A","M1116_1E_1A","C0828_4E_1A","M1088_1E_1A","M1109_1E_1A","M1172_2E_1A","C0915_4E_1A","C0886_4E_1A","C0863_5E_1A","M1344_1E_1A","C1162_3E_1A","C0932_8E_1A","C1068_8E_1A","C1280_4E_1A","C1135_5E_1A","C1271_4E_1A","C1197_2E_1A","M1372_2E_1A","M1344_3E_1A","M1345_3E_1A","C1068_10E_1A","C1155_4E_1A","C1159_4E_1A","M1314_3E_1A","C1273_4E_1A","M1322_3E_1A","C1050_9E_1A","C1220_4E_1A","C1227_3E_1A","C1234_4E_1A","C1339_3E_1A","C1208_4E_1A","C0961_5E_1A","C1007_4E_1A","C1060_4E_1A","M1401_1E_1A","C1116_3E_1A","M1398_1E_1A","C1207_4E_1A","M0514_3E_1A","C1258_4E_1A","C1087_4E_1A","C1237_4E_1A","M1401_2E_1A","C1314_4E_1A","C1326_4E_1A","C1345_4E_1A","C1271_5E_1A","C1050_5E_1A","C1372_3E_1A","C1015_3E_1A","M1273_1E_1A","M1208_3E_1A","M1237_2E_1A","C1087_3E_1A","C0676_3E_1A","M1207_3E_1A","C1061_4E_1A","C0792_4E_1A","M1246_3E_1A","C0774_4E_1A","C0711_6E_1A","M1271_2E_1A","C1135_4E_1A","C0980_5E_1A","M1258_3E_1A","C1177_4E_1A","C1219_3E_1A","C1115_4E_1A","M1280_2E_1A","M1273_3E_1A","C0893_4E_1A","M1339_1E_1A","C1172_4E_1A","C1195_4E_1A","C0723_5E_1A","C0864_6E_1A","C1167_4E_1A","C1099_3E_1A","C0742_5E_1A","M1415_1E_1A","M0891_1E_1A","C0932_3E_1A","M1208_1E_1A","C0794_4E_1A","C1008_4E_1A"])

allmeta = CSV.read("data/metadata/joined.csv", DataFrame)
df = @chain df begin
    leftjoin(select(allmeta, Not([:sample])), on=[:subject, :timepoint])
    filter!(:sample => (s-> in(s, samplessent)), _)
end

boxes = CSV.read("/home/kevin/Desktop/Echo_Metabolomics_Dec20.csv", DataFrame)
boxes = leftjoin(df, select(boxes, Not([:sample, :DOC, :SampleNumber])), on=[:subject, :timepoint])

filter!("metabolomics Box #" => !ismissing, boxes)
sort!(boxes, "metabolomics Box #")
unique!(boxes)

CSV.write("/home/kevin/Desktop/metabolomics_samples.csv", 
   select!(boxes, ["sample", "SampleNumber", "DOC", "subject", "correctedAgeDays", 
                    "Mother_Child", "metabolomics Box #"]))
