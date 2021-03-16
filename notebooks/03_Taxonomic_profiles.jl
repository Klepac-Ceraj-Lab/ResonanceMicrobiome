# # Analysis of Taxonomic Profiles
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

figure1 = Figure(resolution=(1200, 800));

fig1a = figure1[1,1] = Axis(figure1, title="All participants", xlabel=mds_format(all_pco, 1), ylabel=mds_format(all_pco, 2))
scatter!(fig1a, projection(all_pco)[:,1], projection(all_pco)[:,2], 
        color=categorical_colors(all_metadata.ageLabel, ["Prenatal", "1 and under", "1 to 2", "2 and over", missing], colormap[[9, 2, 3, 5, 15]]))
figure1

#- 

kids_species = all_species[:, startswith.(samplenames(all_species), 'C')]
kids_species = kids_species[vec(featuretotals(kids_species) .!= 0), :]

kids_metadata = resonance_metadata(name.(samples(kids_species)))

kids_metadata.ageLabel

kids_pco = pcoa(kids_species)

fig1b = figure1[1,2] = Axis(figure1, title="Children", xlabel=mds_format(kids_pco, 1), ylabel=mds_format(kids_pco, 2))
scatter!(fig1b, projection(kids_pco)[:,1] .* -1, projection(kids_pco)[:,2],
        color=categorical_colors(kids_metadata.ageLabel, ["Prenatal", "1 and under", "1 to 2", "2 and over", missing], colormap[[9, 2, 3, 5, 15]]))

fig1ab_legend = figure1[1,3] = Legend(figure1,
    [
        MarkerElement(color = colormap[9], marker = :circ, strokecolor = :black)
        MarkerElement(color = colormap[2], marker = :circ, strokecolor = :black)
        MarkerElement(color = colormap[3], marker = :circ, strokecolor = :black)
        MarkerElement(color = colormap[5], marker = :circ, strokecolor = :black)
        MarkerElement(color = colormap[15], marker = :circ, strokecolor = :black)
    ],
    ["mom", "1 and under", "1 to 2", "over 2", "missing"])

has_age_species = kids_species[:, map(!ismissing, kids_metadata.correctedAgeDays)]
has_age_metadata = resonance_metadata(name.(samples(has_age_species)))

has_age_metadata.shannon = vec(shannon(has_age_species))

has_age_pco = pcoa(has_age_species)

fig1c = figure1[2, 1:2] = Axis(figure1, xlabel=mds_format(kids_pco, 1), ylabel="Age (years)")

scatter!(fig1c, projection(has_age_pco)[:,1] .* -1, has_age_metadata.correctedAgeDays ./ 365, color=has_age_metadata.shannon, label="Shannon diversity")

fig1c_legend = figure1[2,3] = Colorbar(figure1, halign=:left, limits=extrema(has_age_metadata.shannon), width=25, label="Shannon diversity", )
figure1

#- 
u1_spec = kids_species[:, filter(:correctedAgeDays=> <(365), has_age_metadata).sample]
mid_spec = kids_species[:, filter(:correctedAgeDays=> a-> (365 <= a < 365 * 2), has_age_metadata).sample]
o2_spec = kids_species[:, filter(:correctedAgeDays=> >(365*2), has_age_metadata).sample]

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

function plottopn!(fig, axrow, cp, clust, n, title; kwargs...)
    topnames, topabund = topx(cp, n)
    topabund = topabund[:, clust.order]
    ax = Axis(fig[axrow,1], title=title)
    leg = Legend(fig[axrow, 2],
        [MarkerElement(color = name_dict[topnames[i]], marker = :circle, strokecolor = :black) for i in 1:n+1],
        topnames)
    
    for i in (n+1):-1:1
        @info i, topnames[i]
        v = vec(sum(topabund[1:i, :], dims=1))
        @warn first(v,3)
        barplot!(ax, 1:nsamples(cp), v, color=name_dict[topnames[i]])
    end

    tightlimits!(ax)
    hidexdecorations!(ax)
    fig
end

figure2 = Figure(resolution=(1200,1200));
plottopn!(figure2, 1, u1_spec, u1clust, 10, "Top 10 species, kids under 1 yo")
plottopn!(figure2, 2, mid_spec, midclust, 10, "Top 10 species, kids 1-2 yo")
plottopn!(figure2, 3, o2_spec, o2clust, 10, "Top 10 species, kids over 2 yo")

