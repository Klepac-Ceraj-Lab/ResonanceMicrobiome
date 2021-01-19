using CSV
using DataFrames
using Microbiome
using ECHOAnalysis
using Statistics
using CairoMakie

kneads = String[]
for (root, dirs, files) in walkdir("/Lovelace/echo/analysis/biobakery3/")
    filter!(f-> occursin("kneaddata_read_counts", f), files)
    append!(kneads, joinpath.(root, files))
end
allmeta = CSV.read("data/metadata/joined.csv", DataFrame)

kneaddf = reduce(vcat, CSV.read.(kneads, DataFrame))
kneaddf.final = kneaddf."final pair1" .+ kneaddf."final pair2"
kneaddf.Sample = replace.(replace.(kneaddf.Sample, Ref(r"_S\d+_kneaddata"=>"")), Ref("-"=>"_"))
CSV.write("/home/kevin/Desktop/kneaddata_read_counts.csv", kneaddf)

unirefs = functional_profiles(kind="genefamilies_relab", filefilter=f-> sampleid(stoolsample(basename(f))) in kneaddf.Sample)[1]
unirefs.uniref = replace.(unirefs.func, Ref("UniRef90_"=>""))

milk = CSV.read("data/uniprot/milk.tsv", DataFrame)
milkur = Set(replace.(milk."Entry name", Ref(r"_\w+$"=>"")))

milkunirefs = filter(row-> in(row.uniref, milkur), unirefs)
grps = groupby(unirefs, :sample)

milkdiv = combine(grps, :abundance=>ginisimpson=>:milkgini, :abundance=>shannon=>:milkshannon)
CSV.write("/home/kevin/Desktop/milkdiv.csv", milkdiv)

milkdiv.mc = first.(milkdiv.sample)
grpmc = groupby(milkdiv, :mc)
combine(grpmc, :milkgini=>mean, :milkshannon=>mean)

milkdiv = leftjoin(allmeta, milkdiv, on=:sample)
using Chain

@chain milkdiv begin
    filter(row-> !ismissing(row.correctedAgeDays) && 
                 row.correctedAgeDays < 365 &&
                 !ismissing(row.breastfeeding),_)
    groupby(:breastfeeding)
    combine(:milkshannon=>mean∘skipmissing=>:meanshannon, :milkgini=>mean∘skipmissing=>:meangini)
end

hmo = CSV.read("data/uniprot/uniprot-hmo.tsv", DataFrame)
hmour = Set(replace.(hmo."Entry name", Ref(r"_\w+$"=>"")))
hmo."Entry"

hmounirefs = filter(row-> in(row.uniref, hmour), unirefs)
grps = groupby(hmounirefs, :sample)

hmodiv = combine(grps, :abundance=>ginisimpson=>:hmogini, :abundance=>shannon=>:hmoshannon)
CSV.write("/home/kevin/Desktop/hmodiv.csv", hmodiv)

hmodiv.mc = first.(hmodiv.sample)
grpmc = groupby(hmodiv, :mc)
combine(grpmc, :hmogini=>mean, :hmoshannon=>mean)

hmodiv = leftjoin(allmeta, hmodiv, on=:sample)
using Chain

hmostats = @chain hmodiv begin
    filter(row-> !ismissing(row.correctedAgeDays) && 
                 row.correctedAgeDays < 365 &&
                 !ismissing(row.breastfeeding),_)
    groupby(:breastfeeding)
    combine(:hmoshannon=>mean∘skipmissing=>:meanshannon, :hmogini=>mean∘skipmissing=>:meangini)
end

@chain hmodiv begin
    filter(row-> !ismissing(row.correctedAgeDays) && 
                 row.correctedAgeDays < 365 &&
                 !ismissing(row.breastfeeding),_)
    groupby(:breastfeeding)
end

fig = Figure()

fig_a = fig[1,1] = Axis(fig, title="Exclusive BF")
fig_b = fig[1,2] = Axis(fig, title="Mixed")
fig_c = fig[1,3] = Axis(fig, title="Exclusive Form")

@chain hmodiv begin
    filter(row-> !ismissing(row.correctedAgeDays) && 
                 row.correctedAgeDays < 365 &&
                 !ismissing(row.breastfeeding) &&
                 row.breastfeeding == "exclusive breast", _)
    AbstractPlotting.histogram!(fig_a, collect(skipmissing(_.hmoshannon)))
end

@chain hmodiv begin
    filter(row-> !ismissing(row.correctedAgeDays) && 
                 row.correctedAgeDays < 365 &&
                 !ismissing(row.breastfeeding) &&
                 row.breastfeeding == "mixed", _)
    AbstractPlotting.histogram!(fig_b, collect(skipmissing(_.hmoshannon)))
end

@chain hmodiv begin
    filter(row-> !ismissing(row.correctedAgeDays) && 
                 row.correctedAgeDays < 365 &&
                 !ismissing(row.breastfeeding) &&
                 row.breastfeeding == "exclusive formula", _)
    AbstractPlotting.histogram!(fig_c, collect(skipmissing(_.hmoshannon)))
end

fig