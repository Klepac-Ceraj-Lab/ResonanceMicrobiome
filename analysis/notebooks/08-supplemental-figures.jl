# Supplementary Figures

include("../scripts/startup_loadpackages.jl")
using MakieLayout
using AbstractPlotting
using StatsMakie
using Makie
using ColorSchemes
using FileIO
using StatsBase: midpoints

AbstractPlotting.inline!(false)

@load "analysis/figures/assets/metadata.jld2" allmeta ubothmeta allkidsmeta ukidsmeta allmoms allkids umoms ukids oldkids uboth
allkidsmeta.sample = [String(s) for s in allkidsmeta.sample]
@load "analysis/figures/assets/taxa.jld2" species speciesmds speciesmdsaxes ubothspeciesmds ubothspeciesmdsaxes ukidsspeciesmds ukidsspeciesmdsaxes
@load "analysis/figures/assets/unirefs.jld2" unirefaccessorymds unirefaccessorymdsaxes ubothunirefaccessorymds ubothunirefaccessorymdsaxes ukidsunirefaccessorymds ukidsunirefaccessorymdsaxes
@load "analysis/figures/assets/otherfunctions.jld2" kos kosdiffs kosdm ecs ecsdm pfams pfamsdiffs pfamsdm
@load "analysis/figures/assets/permanovas.jld2" r2 r2m qa allpermanovas species_permanovas unirefaccessory_permanovas kos_permanovas pfams_permanovas
@load "analysis/figures/assets/fsea.jld2" allfsea mdcors
@load "analysis/figures/assets/difs.jld2" speciesdiffs unirefaccessorydiffs kosdiffs pfamsdiffs
@load "analysis/figures/assets/stratkos.jld2" stratkos
@load "analysis/figures/assets/cogquartiles.jld2" quartmeta quartspecies quartspeciesdm quartspeciesmds quartspeciesmdsaxes #quartiletests

allfsea.median = map(median, allfsea.cors)
allmeta.cogAssessment = [x == "None" ? missing : x for x in allmeta.cogAssessment]

set_theme!(
    LAxis = (titlesize=30, xlabelsize=20, ylabelsize=20),
    LLegend = (labelsize=25, markersize=20, patchlabelgap=20)
)

## 

function makelong(diffdict, colors = ColorSchemes.Set3_5.colors, groupsize=4)
    df = DataFrame(xi = Int[], xlabel=String[], y = Float64[], color=[])
    i = 0
    for k1 in sort(keys(diffdict)|> collect)
        subdict = diffdict[k1]
        c = 0
        for k2 in sort(keys(subdict)|> collect)
            i+=1
            c+=1
            dists = subdict[k2]
            append!(df, ((xi = i + floor(Int, (i-0.1)/groupsize), xlabel="$k1 - $k2", y = d, color=colors[c]) for d in dists))
        end
    end
    filter!(row-> row.y != 0, df)
    return df
end
speciesdiffs2 = makelong(speciesdiffs)
unirefaccessorydiffs2 = makelong(unirefaccessorydiffs)
kosdiffs2 = makelong(kosdiffs)
pfamsdiffs2 = makelong(pfamsdiffs)

##

res = (4*300, 3*300)
s1_scene, s1_layout = layoutscene(resolution = res)

axes_diffs = s1_layout[1:2,1:2] = [LAxis(s1_scene, titlevisible=true) for row in 1:2, col in 1:2]

legend_diffs = s1_layout[3,1:2] = LLegend(s1_scene,
    [MarkerElement(marker = :rect, color=ColorSchemes.Set3_5.colors[i], strokecolor=:black) for i in 1:4],
    ["1 and under","1 to 2","2 and over","mom"],
    "Compared to",
    height=Auto(true), tellheight=true, orientation=:horizontal)
s1_scene

s1_layout[1, 1, TopLeft()] = LText(scene, "a", textsize = 40, padding = (0, 0, 10, 0))
s1_layout[1, 2, TopLeft()] = LText(scene, "b", textsize = 40, padding = (0, 0, 10, 0))
s1_layout[2, 1, TopLeft()] = LText(scene, "c", textsize = 40, padding = (0, 0, 10, 0))
s1_layout[2, 2, TopLeft()] = LText(scene, "d", textsize = 40, padding = (0, 0, 10, 0))

##

s1a = let
    plt = nothing
    by(speciesdiffs2, :xi) do df
        plt = boxplot!(axes_diffs[1,1], Data(df),
            :xi, :y,
            color=first(df.color), outliercolor=first(df.color),
            markersize = 7 * AbstractPlotting.px
            )
        end
    plt
end

s1b = let
    plt = nothing
    by(unirefaccessorydiffs2, :xi) do df
        plt = boxplot!(axes_diffs[1,2], Data(df),
            :xi, :y,
            color=first(df.color), outliercolor=first(df.color),
            markersize = 7 * AbstractPlotting.px
            )
        end
    plt
end

s1c = let
    plt = nothing
    by(pfamsdiffs2, :xi) do df
        plt = boxplot!(axes_diffs[2,1], Data(df),
            :xi, :y,
            color=first(df.color), outliercolor=first(df.color),
            markersize = 7 * AbstractPlotting.px
            )
        end
    plt
end

s1d = let
    plt = nothing
    by(kosdiffs2, :xi) do df
        plt = boxplot!(axes_diffs[2,2], Data(df),
            :xi, :y,
            color=first(df.color), outliercolor=first(df.color),
            markersize = 7 * AbstractPlotting.px
            )
        end
    plt
end

axes_diffs[1,1].title = "Species Diffs"
axes_diffs[1,2].title = "Accessory Diffs"
axes_diffs[2,1].title = "Pfams Diffs"
axes_diffs[2,2].title = "KOs Diffs"



for i in 1:4
    axes_diffs[i].xticks = (2.5:5:20, ["1 and under","1 to 2","2 and over","mom"])
    axes_diffs[i].xticksvisible = false
    axes_diffs[i].xgridvisible = false

    ylims!(axes_diffs[i], (0.,1))
end

save("analysis/figures/suppfigure_diffs.jpeg", s1_scene, resolution=res)
s1_scene

##

res=(1200, 900)
s2_scene, s2_layout = layoutscene(resolution = res)
s2_scene

##

s2a = s2_layout[1,1] = LAxis(s2_scene, title="Taxonomic Profiles MDS")

scatter!(s2a, Group(ubothmeta.ageLabel),
        projection(ubothspeciesmds)[:,1], projection(ubothspeciesmds)[:,2],
        grid=false, color=ColorSchemes.Set3_5.colors[1:4],
        markersize = 15 * AbstractPlotting.px,
        strokecolor=:black, strokewidth=1)

s2a.xlabel = "MDS1 ($(round(ubothspeciesmdsaxes[1]*100, digits=2)) %)"
s2a.ylabel = "MDS2 ($(round(ubothspeciesmdsaxes[2]*100, digits=2)) %)"
s2a.xgridvisible = false
s2a.ygridvisible = false
s2a.xticksvisible = false
s2a.yticksvisible = false
s2a.xticklabelsvisible = false
s2a.yticklabelsvisible = false

##

s2b = s2_layout[1,2] = LAxis(s2_scene, title="Uniref90 Accessory MDS")

scatter!(s2b, Group(ubothmeta.ageLabel),
        projection(ubothunirefaccessorymds)[:,1:2],
        grid=false, color=ColorSchemes.Set3_5.colors[1:4],
        markersize = 15 * AbstractPlotting.px,
        strokecolor=:black, strokewidth=1)
s2b.xlabel = "MDS1 ($(round(ubothunirefaccessorymdsaxes[1]*100, digits=2)) %)"
s2b.ylabel = "MDS2 ($(round(ubothunirefaccessorymdsaxes[2]*100, digits=2)) %)"
s2b.xgridvisible = false
s2b.ygridvisible = false
s2b.xticksvisible = false
s2b.yticksvisible = false
s2b.xticklabelsvisible = false
s2b.yticklabelsvisible = false

legend_all_mds = s2_layout[2,1:2] = LLegend(s2_scene,
    [MarkerElement(marker = :circle, color=ColorSchemes.Set3_5.colors[i], strokecolor=:black) for i in 1:4],
    ["1 and under","1 to 2","2 and over","mom"],
    "Subject Age", orientation=:horizontal,
    height=Auto(true), tellheight=true)

##

s2_layout[1, 1, TopLeft()] = LText(s2_scene, "a", textsize = 40, padding = (0, 0, 10, 0))
s2_layout[1, 2, TopLeft()] = LText(s2_scene, "b", textsize = 40, padding = (0, 0, 10, 0))

save("analysis/figures/suppfigure_all_mds.jpeg", s2_scene, resolution=res)
s2_scene

##

res=(800, 800)
s3_scene, s3_layout = layoutscene(resolution = res)
s3_scene

##

s3 = s3_layout[1,1] = LAxis(s3_scene)
pcopri = scatter!(s3, Group(marker=ubothmeta.ageLabel),
        Style(color=log2.(ubothmeta.pcopri .+ (minimum(ubothmeta.pcopri) / 4))), grid=false,
        projection(ubothspeciesmds)[:,1:2],
        markersize = 15 * AbstractPlotting.px, marker=[:utriangle, :hexagon, :rect, :cross],
        strokecolor=:black, strokewidth=1, colormap=:BuPu)

s3.xlabel = "MDS1 ($(round(ubothspeciesmdsaxes[1]*100, digits=2)) %)"
s3.ylabel = "MDS2 ($(round(ubothspeciesmdsaxes[2]*100, digits=2)) %)"
s3.xgridvisible = false
s3.ygridvisible = false
s3.xticksvisible=false
s3.xticklabelsvisible=false
s3.yticksvisible=false
s3.yticklabelsvisible=false

s3_layout[1,2] = LColorbar(s3_scene, pcopri, width=30)
s3_layout[1,2, Left()] = LText(s3_scene, "log2(P. copri abundance)", rotation = pi/2, padding = (0, 5, 0, 0))

marker_legend_pcopri = s3_layout[2,1] = LLegend(s3_scene,
    [
    MarkerElement(marker = :utriangle, color=:white, strokecolor=:black),
    MarkerElement(marker = :hexagon, color=:white, strokecolor=:black),
    MarkerElement(marker = :rect, color=:white, strokecolor=:black),
    MarkerElement(marker = :cross, color=:white, strokecolor=:black)
    ],
    [
    "1 and under",
    "1 to 2",
    "2 and over",
    "mom"
    ],
    "Subject Age", height=Auto(), tellheight=true, tellwidth=false, orientation=:horizontal)
##
save("analysis/figures/suppfigure_pcopri.jpeg", s3_scene, resolution=res)
s3_scene

##

res=(1200, 1200)
s4_scene, s4_layout = layoutscene(resolution = res)
s4_scene
##

cogScore_fsea = s4_layout[1,1] = LAxis(s4_scene, title="Cognitve function score")
plot!(cogScore_fsea, histogram, mdcors[:cogScore])

neocortical_normed_fsea = s4_layout[1,2] = LAxis(s4_scene, title="Neocortical volume")
plot!(neocortical_normed_fsea, histogram, mdcors[:neocortical_normed])

subcortical_normed_fsea = s4_layout[2,1] = LAxis(s4_scene, title="Subcortical volume")
plot!(subcortical_normed_fsea, histogram, mdcors[:subcortical_normed])

limbic_normed_fsea = s4_layout[2,2] = LAxis(s4_scene, title="Limbic volume")
plot!(limbic_normed_fsea, histogram, mdcors[:limbic_normed])

cerebellar_normed_fsea = s4_layout[3,1] = LAxis(s4_scene, title="Cerebellar volume")
plot!(cerebellar_normed_fsea, histogram, mdcors[:cerebellar_normed])



s4_layout[1, 1, TopLeft()] = LText(s4_scene, "a", textsize = 40, padding = (0, 0, 10, 0))
s4_layout[1, 2, TopLeft()] = LText(s4_scene, "b", textsize = 40, padding = (0, 0, 10, 0))
s4_layout[2, 1, TopLeft()] = LText(s4_scene, "c", textsize = 40, padding = (0, 0, 10, 0))
s4_layout[2, 2, TopLeft()] = LText(s4_scene, "d", textsize = 40, padding = (0, 0, 10, 0))
s4_layout[3, 1, TopLeft()] = LText(s4_scene, "b", textsize = 40, padding = (0, 0, 10, 0))

foreach(LAxis, s4_layout) do obj
    tight_ticklabel_spacing!(obj)
    obj.ylabel = "count"
    obj.xlabel = "Pearson correlation"
end

save("analysis/figures/suppfigure_fsea_hist.jpeg", s4_scene, resolution=res)
s4_scene

##
res=(1200, 1600)
s5_scene, s5_layout = layoutscene(resolution=res)
s5_scene

##
brainvols = select(filter(row->!ismissing(row.hires_total), ukidsmeta), 
        [:correctedAgeDays,
        :neocortical_normed, :subcortical_normed, :cerebellar_normed, :limbic_normed,
        :neocortical, :subcortical, :cerebellar, :limbic])
disallowmissing!(brainvols)

neocortical = s5_layout[1,1] = LAxis(s5_scene, title="Neocortex")
scatter!(neocortical, brainvols.correctedAgeDays, brainvols.neocortical, markersize = 10 * AbstractPlotting.px)
neocortical_normed = s5_layout[1,2] = LAxis(s5_scene, title="Neocortex age normalized")
scatter!(neocortical_normed, brainvols.correctedAgeDays, brainvols.neocortical_normed, markersize = 10 * AbstractPlotting.px)

subcortical = s5_layout[2,1] = LAxis(s5_scene, title="Subcortex")
scatter!(subcortical, brainvols.correctedAgeDays, brainvols.subcortical, markersize = 10 * AbstractPlotting.px)
subcortical_normed = s5_layout[2,2] = LAxis(s5_scene, title="Subcortex age normalized")
scatter!(subcortical_normed, brainvols.correctedAgeDays, brainvols.subcortical_normed, markersize = 10 * AbstractPlotting.px)

cerebellar = s5_layout[3,1] = LAxis(s5_scene, title="Cerebellum")
scatter!(cerebellar, brainvols.correctedAgeDays, brainvols.cerebellar, markersize = 10 * AbstractPlotting.px)
cerebellar_normed = s5_layout[3,2] = LAxis(s5_scene, title="Cerebellum age normalized")
scatter!(cerebellar_normed, brainvols.correctedAgeDays, brainvols.cerebellar_normed, markersize = 10 * AbstractPlotting.px)

limbic = s5_layout[4,1] = LAxis(s5_scene, title="Limbic")
scatter!(limbic, brainvols.correctedAgeDays, brainvols.limbic, markersize = 10 * AbstractPlotting.px)
limbic_normed = s5_layout[4,2] = LAxis(s5_scene, title="Limbic age normalized")
scatter!(limbic_normed, brainvols.correctedAgeDays, brainvols.limbic_normed, markersize = 10 * AbstractPlotting.px)
