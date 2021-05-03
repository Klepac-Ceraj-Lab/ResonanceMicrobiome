# # Analysis of MGX Taxonomic Profiles
# 
# The package definied by this repository (`ResonanceMicrobiome`)
# reexports the tools from `Microbiome.jl`, 
# which make it easy to work with microbial taxonomic profiles.

using ResonanceMicrobiome
using Microbiome.MultivariateStats
using CairoMakie
using AbstractPlotting.ColorSchemes
using Clustering

colormap = ColorSchemes.tab20.colors

#-

all_species = taxonomic_profiles(:species)
all_metadata = resonance_metadata(name.(samples(all_species)))

all_pco = pcoa(all_species)

##-

figure1 = Figure(resolution=(800, 800));

fig1a = figure1[1,1] = Axis(figure1, title="All participants", xlabel=mds_format(all_pco, 1), ylabel=mds_format(all_pco, 2))
scatter!(fig1a, projection(all_pco)[:,1], projection(all_pco)[:,2], 
        color=categorical_colors(all_metadata.ageLabel, ["1 and under", "1 to 2", "2 and over", missing], colormap[[2, 3, 5, 15]]))
figure1

#- 

fig1ab_legend = figure1[1,2] = Legend(figure1,
    [
        MarkerElement(color = colormap[2], marker = :circ, strokecolor = :black)
        MarkerElement(color = colormap[3], marker = :circ, strokecolor = :black)
        MarkerElement(color = colormap[5], marker = :circ, strokecolor = :black)
        MarkerElement(color = colormap[15], marker = :circ, strokecolor = :black)
    ],
    ["1 and under", "1 to 2", "over 2", "missing"])

has_age_species = all_species[:, map(!ismissing, all_metadata.correctedAgeDays)]
has_age_metadata = resonance_metadata(name.(samples(has_age_species)))

has_age_metadata.shannon = vec(shannon(has_age_species))

has_age_pco = pcoa(has_age_species)

fig1b = figure1[2, 1] = Axis(figure1, xlabel=mds_format(all_pco, 1), ylabel="Age (years)")

scatter!(fig1b, projection(has_age_pco)[:,1], has_age_metadata.correctedAgeDays ./ 365, color=has_age_metadata.shannon, label="Shannon diversity")

fig1b_legend = figure1[2,2] = Colorbar(figure1, halign=:left, limits=extrema(has_age_metadata.shannon), width=25, label="Shannon diversity", )
CairoMakie.save("figures/03_taxonomic_profiles.svg", figure1)
figure1

#- 

u1_spec = all_species[:, filter(:correctedAgeDays=> <(365), has_age_metadata).sample]
mid_spec = all_species[:, filter(:correctedAgeDays=> a-> (365 <= a < 365 * 2), has_age_metadata).sample]
o2_spec = all_species[:, filter(:correctedAgeDays=> >(365*2), has_age_metadata).sample]

u1dm = braycurtis(u1_spec)
middm = braycurtis(mid_spec)
o2dm = braycurtis(o2_spec)

u1clust = hclust(u1dm, linkage=:complete, branchorder=:optimal)
midclust = hclust(middm, linkage=:complete, branchorder=:optimal)
o2clust = hclust(o2dm, linkage=:complete, branchorder=:optimal)


function topx(cp, n=10)
    totals = vec(featuretotals(cp))
    rows = partialsortperm(totals, 1:n, rev=true)
    top = cp[rows, :]
    other = sum(abundances(cp[collect(1:nfeatures(cp))[Not(rows)], :]), dims=1)

    return (names = vcat(featurenames(top), ["other"]), abundances = vcat(abundances(top), other))
end

allnames = setdiff(union([topx(cp).names for cp in (u1_spec, mid_spec, o2_spec)]...), Set(["other"]))
name_dict = Dict(n => colormap[i] for (i, n) in enumerate(allnames))
name_dict["other"] = ColorSchemes.Greys_3.colors[1]

function plottopn!(fig, axrow, cp, clust, n, title, markerspecies = String[]; kwargs...)
    topnames, topabund = topx(cp, n)
    topabund = topabund[:, clust.order]
    ax = Axis(fig[axrow,1], title=title)
    append!(markerspecies, [topnames[i] for i in 1:n+1])
    
    for i in (n+1):-1:1
        v = vec(sum(topabund[1:i, :], dims=1))
        barplot!(ax, 1:nsamples(cp), v, color=name_dict[topnames[i]])
    end

    tightlimits!(ax)
    hidexdecorations!(ax)
    fig, markerspecies
end

figure2 = Figure(resolution=(1200,1200));
(_, markerspecies) = plottopn!(figure2, 1, u1_spec, u1clust, 10, "Top 10 species, kids under 1 yo")
plottopn!(figure2, 2, mid_spec, midclust, 10, "Top 10 species, kids 1-2 yo", markerspecies)
plottopn!(figure2, 3, o2_spec, o2clust, 10, "Top 10 species, kids over 2 yo", markerspecies)
unique!(markerspecies)
sort!(markerspecies)
fig2_leg = figure2[1:3, 2] = Legend(figure2, [MarkerElement(color = name_dict[m], marker = :rect, strokecolor = :black) for m in markerspecies], [m for m in markerspecies])

CairoMakie.save("figures/03_top_bugs.svg", figure2)
figure2