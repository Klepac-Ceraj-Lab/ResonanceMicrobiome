include("../scripts/startup_loadpackages.jl")


@load "analysis/figures/assets/metadata.jld2" allmeta ubothmeta ukidsmeta allkidsmeta allmoms allkids umoms ukids oldkids uboth
@load "analysis/figures/assets/taxa.jld2" species speciesmds kidsspeciesmds kidsspeciesmdsaxes
@load "analysis/figures/assets/unirefs.jld2" unirefaccessorymds unirefaccessorymdsaxes kidsunirefaccessorymds kidsunirefaccessorymdsaxes
@load "analysis/figures/assets/otherfunctions.jld2" kos kosdiffs kosdm ecs ecsdm pfams pfamsdiffs pfamsdm
@load "analysis/figures/assets/permanovas.jld2" r2 r2m qa allpermanovas species_permanovas unirefaccessory_permanovas kos_permanovas pfams_permanovas
@load "analysis/figures/assets/fsea.jld2" allfsea mdcors
@load "analysis/figures/assets/difs.jld2" speciesdiffs unirefaccessorydiffs kosdiffs pfamsdiffs

allkidsmeta = filter(row-> !ismissing(row.breastfeeding), allkidsmeta)

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

youngkidsmeta = filter(row-> row.correctedAgeDays < 365, allkidsmeta)

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

## 

using Makie
using MakieLayout
using StatsMakie
using ColorSchemes

sig = filter(row-> row.qvalue < 0.1, youngbftests)
for row in eachrow(sig)
    bug = row.species
    youngkidsmeta[!, Symbol(bug)] = collect(vec(occurrences(view(youngkidsspecies, species=[bug]))))
end
youngkidsmeta.x = map(youngkidsmeta.breastfeeding) do bf
    occursin("formula", bf) ? 0 :
    occursin("breast", bf) ? 2 :
    1
end


scene, layout = layoutscene()

bug1 = sig[1, :species]
ery = layout[1,1] = LAxis(scene, title=bug1, titlefont="DejaVu Sans Oblique", titlesize=25)
boxplot!(ery, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug1),
        markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]])
scatter!(ery, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug1),
        markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]], strokewidth=1, strokecolor=:black)

ery.xticks = ([0,1,2], ["Formula", "Mixed", "Breast"])
ery.xlabel = "Feeding"
ery.ylabel = "Relative abundance"

bug2 = sig[2, :species]
cloacae = layout[1,2] = LAxis(scene, title=bug2, titlefont="DejaVu Sans Oblique", titlesize=25)
boxplot!(cloacae, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug2),
        markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]])
scatter!(cloacae, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug2),
        markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]], strokewidth=1, strokecolor=:black)

cloacae.xticks = ([0,1,2], ["Formula", "Mixed", "Breast"])
cloacae.xlabel = "Feeding"
cloacae.ylabel = "Relative abundance"

bug3 = sig[3, :species]
bart = layout[2,1] = LAxis(scene, title=bug3, titlefont="DejaVu Sans Oblique", titlesize=25)
boxplot!(bart, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug3),
        markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]])
scatter!(bart, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug3),
        markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]], strokewidth=1, strokecolor=:black)

bart.xticks = ([0,1,2], ["Formula", "Mixed", "Breast"])
bart.xlabel = "Feeding"
bart.ylabel = "Relative abundance"

bug4 = sig[4, :species]
casse = layout[2,2] = LAxis(scene, title=bug4, titlefont="DejaVu Sans Oblique", titlesize=25)
boxplot!(casse, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug4),
        markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]])
scatter!(casse, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug4),
        markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]], strokewidth=1, strokecolor=:black)

casse.xticks = ([0,1,2], ["Formula", "Mixed", "Breast"])
casse.xlabel = "Feeding"
casse.ylabel = "Relative abundance"


bug5 = sig[5, :species]
torques = layout[3,1] = LAxis(scene, title=bug5, titlefont="DejaVu Sans Oblique", titlesize=25)
boxplot!(torques, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug5),
        markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]])
scatter!(torques, Data(youngkidsmeta), Group(:breastfeeding), :x, Symbol(bug5),
        markersize=AbstractPlotting.px *10, color=ColorSchemes.Accent_5.colors[[3,2,4]], strokewidth=1, strokecolor=:black)

torques.xticks = ([0,1,2], ["Formula", "Mixed", "Breast"])
torques.xlabel = "Feeding"
torques.ylabel = "Relative abundance"

scene