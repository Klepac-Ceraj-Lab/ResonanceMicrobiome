include("../scripts/startup_loadpackages.jl")
using MakieLayout
using AbstractPlotting
using StatsMakie
using Makie
using ColorSchemes
using FileIO
using StatsBase: midpoints
using CairoMakie
CairoMakie.activate!(type="pdf")

AbstractPlotting.inline!(false)

@load "analysis/figures/assets/metadata.jld2" allmeta ubothmeta allkidsmeta ukidsmeta oldkidsmeta allmoms allkids umoms ukids oldkids uboth
allkidsmeta.sample = [String(s) for s in allkidsmeta.sample]
@load "analysis/figures/assets/taxa.jld2" species speciesmds speciesmdsaxes ubothspeciesmds ubothspeciesmdsaxes ukidsspeciesmds ukidsspeciesmdsaxes
@load "analysis/figures/assets/unirefs.jld2" unirefaccessorymds unirefaccessorymdsaxes ubothunirefaccessorymds ubothunirefaccessorymdsaxes ukidsunirefaccessorymds ukidsunirefaccessorymdsaxes
@load "analysis/figures/assets/otherfunctions.jld2" kos kosdiffs kosdm ecs ecsdm pfams pfamsdiffs pfamsdm
@load "analysis/figures/assets/permanovas.jld2" r2 r2m qa allpermanovas species_permanovas unirefaccessory_permanovas kos_permanovas pfams_permanovas
@load "analysis/figures/assets/fsea.jld2" allfsea oldkidsfsea
@load "analysis/figures/assets/difs.jld2" speciesdiffs unirefaccessorydiffs kosdiffs pfamsdiffs
@load "analysis/figures/assets/stratkos.jld2" stratkos
@load "analysis/figures/assets/cogquartiles.jld2" quartmeta quartspecies quartspeciesdm quartspeciesmds quartspeciesmdsaxes #quartiletests

##

allfsea.median = map(median, allfsea.cors)
oldkidsfsea.median = map(median, oldkidsfsea.cors)
allmeta.cogAssessment = [(ismissing(x) || x == "None") ? missing : x for x in allmeta.cogAssessment]

set_theme!(
    LAxis = (titlesize=30, xlabelsize=20, ylabelsize=20),
    LLegend = (labelsize=25, markersize=20, patchlabelgap=20)
)

##

ebf = findall(x-> !ismissing(x) && x == "exclussive breast", allkidsmeta.breastfeeding)
eff = findall(x-> !ismissing(x) && x == "exclussive formula", allkidsmeta.breastfeeding)

allkidsspecies = view(species, sites=allkidsmeta.sample)
@assert sitenames(allkidsspecies) == allkidsmeta.sample

bftests = DataFrame()

for (i, row) in enumerate(eachrow(occurrences(allkidsspecies)))
    row = vec(row)
    sp = featurenames(allkidsspecies)[i]
    bf = row[ebf]
    ff = row[eff]
    mwu = MannWhitneyUTest(bf, ff)
    assoc = mean(bf) > mean(ff) ? "-" : "+"
    push!(bftests, (
        species=sp,
        association=assoc,
        median_lower = median(bf),
        median_upper = median(ff),
        nsamples = count(>(0), row),
        pvalue = pvalue(mwu)
    ))
end
bftests[!,:qvalue] = adjust(bftests.pvalue, BenjaminiHochberg())
sort!(bftests, :qvalue)
CSV.write("analysis/bftests.csv", bftests)

##

youngkidsmeta = filter(row-> row.correctedAgeDays < (365 / 2), allkidsmeta)

ebf = findall(x-> !ismissing(x) && x == "exclussive breast", youngkidsmeta.breastfeeding)
eff = findall(x-> !ismissing(x) && x == "exclussive formula", youngkidsmeta.breastfeeding)

youngkidsspecies = view(species, sites=youngkidsmeta.sample)
@assert sitenames(youngkidsspecies) == youngkidsmeta.sample

youngbftests = DataFrame()

for (i, row) in enumerate(eachrow(occurrences(youngkidsspecies)))
    row = vec(row)
    sp = featurenames(youngkidsspecies)[i]
    bf = row[ebf]
    ff = row[eff]
    mwu = MannWhitneyUTest(bf, ff)
    assoc = mean(bf) > mean(ff) ? "-" : "+"
    push!(youngbftests, (
        species=sp,
        association=assoc,
        median_lower = median(bf),
        median_upper = median(ff),
        nsamples = count(>(0), row),
        pvalue = pvalue(mwu)
    ))
end
youngbftests[!,:qvalue] = adjust(youngbftests.pvalue, BenjaminiHochberg())
sort!(youngbftests, :qvalue)
CSV.write("analysis/youngbftests.csv", youngbftests)

using BiobakeryUtils
using Distances

youngkidsspeciesdm = pairwise(BrayCurtis(), youngkidsspecies)
youngbfperm = permanova(youngkidsspeciesdm, [ismissing(x) ? missing : string(x) for x in youngkidsmeta.breastfeeding], label="species")


## 

using Makie
using MakieLayout
using StatsMakie
using ColorSchemes

# sig = filter(row-> row.qvalue < 0.2, youngbftests)
# for row in eachrow(sig)
#     bug = row.species
#     youngkidsmeta[!, Symbol(bug)] = collect(vec(occurrences(view(youngkidsspecies, species=[bug]))))
# end
# youngkidsmeta.x = map(youngkidsmeta.breastfeeding) do bf
#     ismissing(bf) ? -1 : 
#     occursin("formula", bf) ? 0 :
#     occursin("breast", bf) ? 2 :
#     1
# end


# scene, layout = layoutscene()

# bug1 = sig[1, :species]
# ery = layout[1,1] = LAxis(scene, title=bug1, titlefont="DejaVu Sans Oblique", titlesize=25)
# boxplot!(ery, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug1),
#         markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]])
# scatter!(ery, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug1),
#         markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]], strokewidth=1, strokecolor=:black)

# ery.xticks = ([0,1,2], ["Formula", "Mixed", "Breast"])
# ery.xlabel = "Feeding"
# ery.ylabel = "Relative abundance"

# bug2 = sig[2, :species]
# cloacae = layout[1,2] = LAxis(scene, title=bug2, titlefont="DejaVu Sans Oblique", titlesize=25)
# boxplot!(cloacae, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug2),
#         markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]])
# scatter!(cloacae, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug2),
#         markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]], strokewidth=1, strokecolor=:black)

# cloacae.xticks = ([0,1,2], ["Formula", "Mixed", "Breast"])
# cloacae.xlabel = "Feeding"
# cloacae.ylabel = "Relative abundance"

# bug3 = sig[3, :species]
# bart = layout[2,1] = LAxis(scene, title=bug3, titlefont="DejaVu Sans Oblique", titlesize=25)
# boxplot!(bart, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug3),
#         markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]])
# scatter!(bart, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug3),
#         markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]], strokewidth=1, strokecolor=:black)

# bart.xticks = ([0,1,2], ["Formula", "Mixed", "Breast"])
# bart.xlabel = "Feeding"
# bart.ylabel = "Relative abundance"

# bug4 = sig[4, :species]
# casse = layout[2,2] = LAxis(scene, title=bug4, titlefont="DejaVu Sans Oblique", titlesize=25)
# boxplot!(casse, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug4),
#         markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]])
# scatter!(casse, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug4),
#         markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]], strokewidth=1, strokecolor=:black)

# casse.xticks = ([0,1,2], ["Formula", "Mixed", "Breast"])
# casse.xlabel = "Feeding"
# casse.ylabel = "Relative abundance"


# bug5 = sig[5, :species]
# torques = layout[3,1] = LAxis(scene, title=bug5, titlefont="DejaVu Sans Oblique", titlesize=25)
# boxplot!(torques, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug5),
#         markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]])
# scatter!(torques, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug5),
#         markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]], strokewidth=1, strokecolor=:black)

# torques.xticks = ([0,1,2], ["Formula", "Mixed", "Breast"])
# torques.xlabel = "Feeding"
# torques.ylabel = "Relative abundance"

# scene

## < 1 yo
using StatsBase
using Distances
using MultivariateStats

AbstractPlotting.inline!(false)

under1meta = filter(row-> !ismissing(row.correctedAgeDays) && row.correctedAgeDays < 365, allkidsmeta)
under1meta.bf = categorical(map(x -> ismissing(x) ? "unknown" :
                            x == "mixed" ? "mixed" :
                            x == "exclussive formula" ? "formula" :
                            x == "exclussive breast" ? "breast" :
                            error("unknown value $x"),
                        under1meta.breastfeeding))

bftypes = ["unknown", "mixed", "formula", "breast"]
levels!(under1meta.bf, bftypes)
under1ns = countmap(under1meta.bf)

under1species = view(species, sites=under1meta.sample)
under1ecs = view(ecs, sites=under1meta.sample)
bfcolors = ColorSchemes.Set3_6.colors[3:end]
under1meta.bif = vec(sum(occurrences(view(under1species, species=map(x-> occursin("Bifidobacterium",x), speciesnames(under1species)))), dims=1))
under1meta.glyc = vec(sum(occurrences(view(under1ecs, species=map(x-> occursin(r"3\.2\.1\.",x), speciesnames(under1ecs)))), dims=1))


under1dm = pairwise(BrayCurtis(),under1species)
under1pco = fit(MDS, under1dm, distances=true)
under1pcoax = [v / sum(eigvals(under1pco)) for v in eigvals(under1pco)]

##

scene, layout = layoutscene(resolution=(1200,900))

bifido = layout[1,1] = LAxis(scene)
boxplot!(bifido, Data(under1meta), levelcode.(under1meta.bf), :bif, markersize = 15 * AbstractPlotting.px, color=bfcolors)
bifido.xticks = (1:4, ["$t\n(n=$(under1ns[t]))" for t in bftypes])
bifido.ylabel = "Relative abundance of Bifidobacterium"

save("/Users/ksb/Desktop/plotting.pdf", scene)


under1pcoa_ax = layout[1,2] = LAxis(scene)
under1grouppcoa = scatter!(under1pcoa_ax, Group(under1meta.bf), projection(under1pco)[:,1:2], 
                markersize = 15 * AbstractPlotting.px,
                color=ColorSchemes.Set3_6.colors[[3,4,5,6]])

under1pcoa_ax.xlabel = "MDS1 ($(round(under1pcoax[1]*100, digits=2)) %)"
under1pcoa_ax.ylabel = "MDS2 ($(round(under1pcoax[2]*100, digits=2)) %)"

bfmarkers = [MarkerElement(marker = :circle, color=bfcolors[i], strokecolor=:black) for i in 1:4]
under1pcoa_leg = layout[1,3] = LLegend(scene, bfmarkers, bftypes)

glyc = layout[2,1] = LAxis(scene)
boxplot!(glyc, Data(under1meta), levelcode.(under1meta.bf), :glyc, markersize = 15 * AbstractPlotting.px, color=bfcolors)
glyc.xticks = (1:4, ["$t\n(n=$(under1ns[t]))" for t in bftypes])
glyc.ylabel = "Relative abundance of Glycosidase"



under1agepco_ax = layout[2,2] = LAxis(scene)
under1agepcoa = scatter!(under1agepco_ax, Style(color=Float64.(under1meta.correctedAgeDays)), projection(under1pco)[:,1:2], 
                markersize = 15 * AbstractPlotting.px)
under1agepco_ax.xlabel = "MDS1 ($(round(under1pcoax[1]*100, digits=2)) %)"
under1agepco_ax.ylabel = "MDS2 ($(round(under1pcoax[2]*100, digits=2)) %)"
under1agepco_legend = layout[2,3] = LColorbar(scene, under1agepcoa, width=30)
layout[2, 3, Left()] = LText(scene, "Age (days)", textsize = 25, rotation = pi/2, padding = (0, 5, 0, 0))


layout[0,1:3] = LText(scene, "Kids < 1yo, n=$(nrow(under1meta))", textsize=30)
save("/Users/ksb/Desktop/under1.pdf", scene)
scene

##

## < 6 mo
using StatsBase
using Distances
using MultivariateStats

under6mometa = filter(row-> !ismissing(row.correctedAgeDays) && row.correctedAgeDays < (365/2) && !ismissing(row.breastfeeding), allkidsmeta)
under6mometa.bf = categorical(map(x -> ismissing(x) ? "unknown" :
                            x == "mixed" ? "mixed" :
                            x == "exclussive formula" ? "formula" :
                            x == "exclussive breast" ? "breast" :
                            error("unknown value $x"),
                        under6mometa.breastfeeding))

bftypes = ["mixed", "formula", "breast"]
levels!(under6mometa.bf, bftypes)
under6mons = countmap(under6mometa.bf)

under6mospecies = view(species, sites=under6mometa.sample)
under6moecs = view(ecs, sites=under6mometa.sample)
bfcolors = ColorSchemes.Set3_6.colors[4:end]
under6mometa.bif = vec(sum(occurrences(view(under6mospecies, species=map(x-> occursin("Bifidobacterium",x), speciesnames(under6mospecies)))), dims=1))
under6mometa.glyc = vec(sum(occurrences(view(under6moecs, species=map(x-> occursin(r"3\.2\.1\.",x), speciesnames(under6moecs)))), dims=1))


under6modm = pairwise(BrayCurtis(),under6mospecies)
under6mopco = fit(MDS, under6modm, distances=true)
under6mopcoax = [v / sum(eigvals(under6mopco)) for v in eigvals(under6mopco)]

##

scene, layout = layoutscene(resolution=(1200,900))

bifido = layout[1,1] = LAxis(scene)
boxplot!(bifido, Data(under6mometa), levelcode.(under6mometa.bf), :bif, markersize = 15 * AbstractPlotting.px, color=bfcolors)
bifido.xticks = (1:3, ["$t\n(n=$(under6mons[t]))" for t in bftypes])
bifido.ylabel = "Relative abundance of Bifidobacterium"

save("/Users/ksb/Desktop/plotting.pdf", scene)

##
under6mopcoa_ax = layout[1,2] = LAxis(scene)
groupunder6mopcoa = scatter!(under6mopcoa_ax, Group(under6mometa.bf), projection(under6mopco)[:,1:2], 
                markersize = 15 * AbstractPlotting.px,
                color=ColorSchemes.Set3_6.colors[[4,5,6]])

under6mopcoa_ax.xlabel = "MDS1 ($(round(under6mopcoax[1]*100, digits=2)) %)"
under6mopcoa_ax.ylabel = "MDS2 ($(round(under6mopcoax[2]*100, digits=2)) %)"

bfmarkers = [MarkerElement(marker = :circle, color=bfcolors[i], strokecolor=:black) for i in 1:3]
pcoa_leg = layout[1,3] = LLegend(scene, bfmarkers, bftypes)

##
# glyc = layout[2,1] = LAxis(scene)
# boxplot!(glyc, Data(under6mometa), levelcode.(under6mometa.bf), :glyc, markersize = 15 * AbstractPlotting.px, color=bfcolors)
# glyc.xticks = (1:3, ["$t\n(n=$(under6mons[t]))" for t in bftypes])
# glyc.ylabel = "Relative abundance of Glycosidase"

##
bcd = layout[2,1] = LAxis(scene)

under6mobfidx = findall(f-> f == "breast", under6mometa.bf)
under6moffidx = findall(f-> f == "formula", under6mometa.bf)

u6mobfdist = filter(!=(0.), vec(under6modm[under6mobfidx,under6mobfidx]))
u6moffdist = vec(under6modm[under6moffidx,under6moffidx])
u6mobfffdist = vec(under6modm[under6mobfidx,under6moffidx])

plotxvals = [fill(1, length(u6mobfdist)); fill(2, length(u6moffdist)); fill(3, length(u6mobfffdist))]
plotyvals = [u6mobfdist; u6moffdist; u6mobfffdist]

boxplot!(bcd, plotxvals, plotyvals, markersize = 15 * AbstractPlotting.px, color=ColorSchemes.seaborn_bright6[1:3])
bcd.xticks = (1:3, ["BF:BF", "FF:FF", "BF:FF"])
bcd.ylabel = "Bray Curtis Dissimilarity"

##
bfcolors
under6moagepco_ax = layout[2,2] = LAxis(scene)
under6moagepcoa = scatter!(under6moagepco_ax, Style(color=Float64.(under6mometa.correctedAgeDays)), projection(under6mopco)[:,1:2], 
                markersize = 15 * AbstractPlotting.px)
under6moagepco_ax.xlabel = "MDS1 ($(round(under6mopcoax[1]*100, digits=2)) %)"
under6moagepco_ax.ylabel = "MDS2 ($(round(under6mopcoax[2]*100, digits=2)) %)"
under6moagepco_legend = layout[2,3] = LColorbar(scene, under6moagepcoa, width=30)
layout[2, 3, Left()] = LText(scene, "Age (days)", textsize = 25, rotation = pi/2, padding = (0, 5, 0, 0))


layout[0,1:3] = LText(scene, "Kids < 6mo, n=$(nrow(under6mometa))", textsize=30)
save("/Users/ksb/Desktop/under6mo.pdf", scene)
scene

##

top = filterabund(under1species, 10)
topspec = featurenames(top)

for sp in topspec
    under1meta[!, Symbol(sp)] = collect(vec(occurrences(view(top, species=[sp]))))
end

under1hcl = hclust(under1dm, linkage=:average, branchorder=:optimal)
long = stack(under1meta, Symbol.(topspec), [:sample, :bf])
long.sample = categorical(long.sample)
levels!(long.sample, String.(under1meta.sample)[under1hcl.order])

##

scene, layout = layoutscene()
topab = layout[1,1] = LAxis(scene)
unique(long.variable) |> length
barplot!(topab, Position.stack, Data(long), Group(color=:variable), :sample, :value, 
    color=ColorSchemes.Set3_11.colors)
tightlimits!(topab)

spmarkers = [MarkerElement(marker = :rect, color=ColorSchemes.Set3_11.colors[i], strokecolor=:black) for i in 1:11]
leg = layout[1,2] = LLegend(scene, spmarkers, topspec)

hm = layout[2,1] = LAxis(scene, height=40)
hmplot = heatmap!(hm, hcat([under1meta.correctedAgeDays[i] for i in under1hcl.order]), interpolate=false)
tightlimits!(hm)
hidedecorations!(hm)

cb = layout[2,2] = LColorbar(scene, hmplot, vertical=false, label="Age (days)")
layout[0,1:2] = LText(scene, "Kids < 1yo, n=$(nrow(under6mometa))", textsize=30)
save("/Users/ksb/Desktop/under1topspec.pdf", scene)
scene

##

top = filterabund(under6mospecies, 15)
topspec = featurenames(top)

for sp in topspec
    under6mometa[!, Symbol(sp)] = collect(vec(occurrences(view(top, species=[sp]))))
end

under6mohcl = hclust(under6modm, linkage=:average, branchorder=:optimal)
long = stack(under6mometa, Symbol.(topspec), [:sample, :bf, :correctedAgeDays])
levels(long.bf)
sort!(long, [:bf, :correctedAgeDays, :sample])
unique(long.sample)

long.sample = categorical(long.sample)
levels!(long.sample, unique(long.sample))
levelcode.(long.sample)
under6mometa.sample = categorical(under6mometa.sample)
levels!(under6mometa.sample, unique(long.sample))

srt = sortperm(under6mometa.sample)
##

scene, layout = layoutscene()
topab = layout[1,1] = LAxis(scene)

unique(long.variable) |> length
barplot!(topab, Position.stack, Data(long), Group(color=:variable), :sample, :value,
    color=ColorSchemes.tab20.colors[1:16])
tightlimits!(topab)

spmarkers = [MarkerElement(marker = :rect, color=ColorSchemes.tab20.colors[i], strokecolor=:black) for i in 1:16]
leg = layout[1,2] = LLegend(scene, spmarkers, topspec)

agehm = layout[2,1] = LAxis(scene, height=40)

hmplot = heatmap!(agehm, reshape(under6mometa[srt,:correctedAgeDays], length(srt), 1), interpolate=false)
tightlimits!(agehm)
hidedecorations!(agehm)

bfann = layout[3,1] = LAxis(scene, height=40)
for (x, i) in enumerate(srt)
    c = under6mometa.bf[i] == "breast" ? bfcolors[3] :
        under6mometa.bf[i] == "formula" ? bfcolors[2] :
        bfcolors[1]

    (start, stop) = (x, x+1)
    poly!(bfann, Point2f0[(start,0), (stop, 0), (stop, 1), (start, 1)], color = c)
end

tightlimits!(bfann)
hidedecorations!(bfann)

bfleg = layout[3,2] = LLegend(scene, bfmarkers, bftypes, orientation=:horizontal)

cb = layout[2,2] = LColorbar(scene, hmplot, vertical=false, label="Age (days)")
layout[0,1:2] = LText(scene, "Kids < 6mo, n=$(nrow(under6mometa))", textsize=30)
save("/Users/ksb/Desktop/under6motopspec.pdf", scene)
scene


## Demographics
levels(under6mometa.simple_race)
n = nrow(under6mometa)
mixedrace = occursin.(r"[Mm]ixed", under6mometa.simple_race) |> sum
whiterace = occursin.(r"[Ww]hite", under6mometa.simple_race) |> sum
africanrace = occursin.(r"African", under6mometa.simple_race) |> sum
unknownrace = occursin.(r"(Unknown|Decline)", under6mometa.simple_race) |> sum

println("Percent non-white: ", round((n - whiterace) / (n - unknownrace) * 100, digits=1))
println("Percent mixed: ", round(mixedrace / (n - unknownrace) * 100, digits=1))
println("Percent African American: ", round(africanrace / (n - unknownrace) * 100, digits=1))
describe(under6mometa.correctedAgeDays)
describe(under6mometa.mother_HHS)
describe(under6mometa.simple_race)

usixallmeta = filter(row-> !ismissing(row.correctedAgeDays) && row.correctedAgeDays < (365/2), allkidsmeta)

n = nrow(usixallmeta)
mixedrace = occursin.(r"[Mm]ixed", usixallmeta.simple_race) |> sum
whiterace = occursin.(r"[Ww]hite", usixallmeta.simple_race) |> sum
africanrace = occursin.(r"African", usixallmeta.simple_race) |> sum
unknownrace = occursin.(r"(Unknown|Decline)", usixallmeta.simple_race) |> sum

println("Percent non-white: ", round((n - whiterace) / (n - unknownrace) * 100, digits=1))
println("Percent mixed: ", round(mixedrace / (n - unknownrace) * 100, digits=1))
println("Percent African American: ", round(africanrace / (n - unknownrace) * 100, digits=1))
print("correctedAgeDays "); describe(usixallmeta.correctedAgeDays)
print("mother_HHS "); describe(usixallmeta.mother_HHS)
describe(usixallmeta.simple_race)

## HMO taxa

usixallspecies = view(species, sites=usixallmeta.sample) |> copy

bugs = [
    "Bifidobacterium_bifidum",
    "Bifidobacterium_breve",
    "Bifidobacterium_longum",
    "Bifidobacterium_kashiwanohense",
    "Bifidobacterium_pseudocatenulatum"
]

open("hmobugs.csv","w")  do io 
    specs = speciesnames(usixallspecies)
    println(io, "bug,min,max,mean,nonzero-min,nonzero-max,nonzero-mean")

    for bug in bugs
        if !in(bug, specs)
            @warn "$bug not found in samples"
            continue
        end
        abs = vec(occurrences(view(usixallspecies, species=[bug])))
        nz = filter(!=(0), abs)
        println(io, join([bug,minimum(abs),maximum(abs),mean(abs),minimum(nz),maximum(nz),mean(nz)], ","))
    end
end

lacto = view(usixallspecies, species=map(s-> occursin("Lacto", s), speciesnames(usixallspecies)))
sum(occurrences(lacto))
for bug in speciesnames(lacto)
    ab = collect(vec(occurrences(view(lacto, species=[bug]))))
    sum(ab) > 0 || continue
    usixallmeta[!, Symbol(bug)] = ab
end
lactometa = select(usixallmeta, [:subject, :timepoint, :sample, names(usixallmeta)[end-nspecies(lacto):end]...])

for n in names(lactometa)[4:end]
    print(n, " ")
    describe(lactometa[!, n])
end

sum(occurrences(usixallspecies)[:,1])