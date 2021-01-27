# # Analysis of Taxonomic Profiles
# 
# The package definied by this repository (`ResonanceMicrobiome`)
# reexports the tools from `Microbiome.jl`, 
# which make it easy to work with microbial taxonomic profiles.

using ResonanceMicrobiome
using Microbiome.MultivariateStats
using CairoMakie
using AbstractPlotting.ColorSchemes

function labelage(df)
    map(eachrow(df)) do row
        if startswith(row.sample, "M")
            "mom"
        elseif ismissing(row.correctedAgeDays)
            missing
        elseif row.correctedAgeDays < 365
            "1 and under"
        elseif row.correctedAgeDays < 2 * 365
            "1 to 2"
        else
            "over 2"
        end
    end
end

colormap = ColorSchemes.tab20.colors

#-

all_species = taxonomic_profiles(:species)
all_metadata = resonance_metadata(name.(samples(all_species)))
all_metadata.ageLabel = labelage(all_metadata)

all_pco = pcoa(all_species)

##-

figure1 = Figure(resolution=(1200, 800));

fig1a = figure1[1,1] = Axis(figure1, title="All participants", xlabel=mds_format(all_pco, 1), ylabel=mds_format(all_pco, 2))
scatter!(fig1a, projection(all_pco)[:,1], projection(all_pco)[:,2], 
        color=categorical_colors(all_metadata.ageLabel, ["mom", "1 and under", "1 to 2", "over 2", missing], colormap[[9, 2, 3, 5, 15]]))
figure1

#- 

kids_species = all_species[:, startswith.(samplenames(all_species), 'C')]
kids_species = kids_species[vec(featuretotals(kids_species) .!= 0), :]

kids_metadata = resonance_metadata(name.(samples(kids_species)))
kids_metadata.ageLabel = labelage(kids_metadata)

kids_pco = pcoa(kids_species)

fig1b = figure1[1,2] = Axis(figure1, title="Children", xlabel=mds_format(kids_pco, 1), ylabel=mds_format(kids_pco, 2))
scatter!(fig1b, projection(kids_pco)[:,1] .* -1, projection(kids_pco)[:,2],
        color=categorical_colors(kids_metadata.ageLabel, ["mom", "1 and under", "1 to 2", "over 2", missing], colormap[[9, 2, 3, 5, 15]]))
figure1

fig1ab_legend = figure1[1,3] = Legend(figure1,
    [
        MarkerElement(color = colormap[9], marker = :circ, strokecolor = :black)
        MarkerElement(color = colormap[2], marker = :circ, strokecolor = :black)
        MarkerElement(color = colormap[3], marker = :circ, strokecolor = :black)
        MarkerElement(color = colormap[5], marker = :circ, strokecolor = :black)
        MarkerElement(color = colormap[15], marker = :circ, strokecolor = :black)
    ],
    ["mom", "1 and under", "1 to 2", "over 2", "missing"])
figure1

has_age_species = kids_species[:, map(!ismissing, kids_metadata.correctedAgeDays)]
has_age_metadata = resonance_metadata(name.(samples(has_age_species)))

has_age_metadata.ageLabel = labelage(has_age_metadata)
has_age_metadata.shannon = vec(shannon(has_age_species))

has_age_pco = pcoa(has_age_species)

fig1c = figure1[2, 1:2] = Axis(figure1, xlabel=mds_format(kids_pco, 1), ylabel="Age (years)")

scatter!(fig1c, projection(has_age_pco)[:,1] .* -1, has_age_metadata.correctedAgeDays ./ 365, color=has_age_metadata.shannon, label="Shannon diversity")

fig1c_legend = figure1[2,3] = Colorbar(figure1, halign=:left, limits=extrema(has_age_metadata.shannon), width=25, label="Shannon diversity", )
figure1