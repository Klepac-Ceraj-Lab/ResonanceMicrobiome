using ResonanceMicrobiome
using BiobakeryUtils
using CairoMakie
# using AlgebraOfGraphics
using AbstractPlotting.ColorSchemes
using Statistics

colormap = ColorSchemes.tab20.colors

#-

all_metadata = resonance_metadata()
all_metadata.correctedAgeMonths = all_metadata.correctedAgeDays ./ 365 .* 12
kids_metadata = filter(row-> !ismissing(row.correctedAgeDays), all_metadata)
has_bf = filter(row-> !ismissing(row.breastfeeding), all_metadata)

figure1 = Figure(resolution=(1200,1600))

fig1a = Axis(figure1[1,1], title="Subject types", ylabel = "Count")
fig1b = Axis(figure1[1,2], title="Child ages", ylabel = "Count")
fig1c = Axis(figure1[2,1], ylabel = "Count")
fig1d = Axis(figure1[2,2], xlabel = "Age (months)", ylabel = "Count")
fig1c_leg = Legend(figure1[2,3], [
    MarkerElement(color = colormap[1], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[7], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[9], marker = :circle, strokecolor = :black),
    ], ["Formula", "Breastmilk", "Mixed"], "Diet")

fig1e = Axis(figure1[3,1], title = "Race / ethnicity", ylabel = "Count")
fig1f = Axis(figure1[3,2], title = "SES", ylabel = "Count")
fig1g = Axis(figure1[4,1], title = "WHO BMI", ylabel = "BMI (z-score)", xlabel="Age (months)")
fig1h = Axis(figure1[4,2], title = "Cognitive score", ylabel = "Score", xlabel="Age (months)")
fig1h_leg = Legend(figure1[4,3], [
    MarkerElement(color = colormap[10], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[13], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[17], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[19], marker = :circle, strokecolor = :black),
    ], ["Mullen", "Bayleys", "WPPSI", "WISC"], "Assessment")

fig1i = Axis(figure1[5,2], ylabel = "Count")
fig1i_leg = Legend(figure1[5,3], [
    MarkerElement(color = colormap[9], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[14], marker = :circle, strokecolor = :black),
    ], ["Cesarean", "Vaginal"], "Delivery method")
    

barplot!(fig1a, [1,2], [count(x-> x === ("Prenatal"), all_metadata.ageLabel), count(x-> x !==("Prenatal"), all_metadata.ageLabel)], color=:gray)
fig1a.xticks = ([1,2], ["moms", "kids"])

hist!(fig1b, kids_metadata.correctedAgeMonths, color=:gray)

barplot!(fig1c, [1,2,3], [
                count(x-> x === ("exclusive formula"), all_metadata.breastfeeding),
                count(x-> x === ("exclusive breast"), all_metadata.breastfeeding),
                count(x-> x === ("mixed"), all_metadata.breastfeeding)],
                color = colormap[[1,7,9]])
fig1c.xticks = ([1,2,3], ["Formula", "Breastmilk", "Mixed"])
fig1c.xticklabelsize = 20

hist!(fig1d, filter(row-> row.breastfeeding == "exclusive formula", has_bf).correctedAgeMonths, color=colormap[1])
hist!(fig1d, filter(row-> row.breastfeeding == "exclusive breast", has_bf).correctedAgeMonths, color=colormap[7])
hist!(fig1d, filter(row-> row.breastfeeding == "mixed", has_bf).correctedAgeMonths, color=colormap[9])

barplot!(fig1e, 1:6,
    [count(r-> !ismissing(r) && occursin("White", r), kids_metadata.simple_race),
    count(r-> !ismissing(r) && occursin("Black", r), kids_metadata.simple_race),
    count(r-> !ismissing(r) && occursin("Asian", r), kids_metadata.simple_race),
    count(r-> !ismissing(r) && occursin("Native", r), kids_metadata.simple_race),
    count(r-> !ismissing(r) && occursin("Mixed", r), kids_metadata.simple_race),
    count(r-> ismissing(r) || occursin("Unknown", r) || occursin("Decline", r), kids_metadata.simple_race)],
    )
fig1e.xticks = (1:6, ["White", "Black", "Asian", "Native", "Mixed", "Unknown"])
fig1e.xticklabelrotation = pi/4

barplot!(fig1f, 3:7, map(n-> count(==(n), skipmissing(kids_metadata.mother_HHS)), 3:7))
fig1f.xticks = 3:7

hasbmi = .!ismissing.(kids_metadata.WHO_zbmi)
scatter!(fig1g, collect(skipmissing(kids_metadata[hasbmi, :correctedAgeMonths])), collect(skipmissing(kids_metadata[hasbmi, :WHO_zbmi])))

hascog = .!ismissing.(kids_metadata.cogScore)
scatter!(fig1h, collect(skipmissing(kids_metadata[hascog, :correctedAgeMonths])), collect(skipmissing(kids_metadata[hascog, :cogScore])), 
        color=categorical_colors(kids_metadata[hascog, :cogAssessment], ["Mullen", "Bayleys", "WPPSI", "WISC"], colormap[[10,13,17,19]]))

barplot!(fig1i, [1,2], [
    count(x-> x === ("Cesarean"), unique(all_metadata, :subject).birthType),
    count(x-> x === ("Vaginal"), unique(all_metadata, :subject).birthType)],
    color = colormap[[9,14]])
fig1i.xticks = ([1,2], ["Cesarean", "Vaginal"])
fig1i.xticklabelsize = 20

figure1

CairoMakie.save("figures/05_data_summaries.svg", figure1)

#- 

has_16S = filter(row-> !ismissing(row."16S_batch"), kids_metadata)
figure2 = Figure(resolution=(1200,800))

age16S = Axis(figure2[1,1], title = "age distribution for 16S")
hist!(age16S, has_16S.correctedAgeMonths)
figure2

describe(has_16S.correctedAgeMonths)
count(a -> 4 < a < 6, has_16S.correctedAgeMonths)
count(a -> 11 < a < 13, has_16S.correctedAgeMonths)
filter(row -> 4 < row.correctedAgeMonths < 6, has_16S).subject ∩ filter(row -> 11 < row.correctedAgeMonths < 13, has_16S).subject

size(all_metadata)
size(kids_metadata)

#- 

figure3a = Figure(resolution=(800,600))
figure3b = Figure(resolution=(800,600))

(lq, uq) = quantile(kids_metadata[hascog, :cogScore], [0.25, 0.75])
kids_metadata.cogQuant = map(c-> ismissing(c)    ? missing : 
                                    c < lq       ? 1       :
                                    lq <= c < uq ? 2 : 3, kids_metadata.cogScore)

fig3a = Axis(figure3a[1,1], title="Cognitive score by test", xlabel="Age (years)")
leg3a = Legend(figure3a[1,2], [
    MarkerElement(color = colormap[10], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[13], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[17], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[19], marker = :circle, strokecolor = :black),
    ], ["Mullen", "Bayleys", "WPPSI", "WISC"], "Assessment", width = 150)
plot!(fig3a,
        collect(skipmissing(kids_metadata[hascog, :correctedAgeMonths])) ./ 12,
        collect(skipmissing(kids_metadata[hascog, :cogScore])), 
        color=categorical_colors(kids_metadata[hascog, :cogAssessment], ["Mullen", "Bayleys", "WPPSI", "WISC"], colormap[[10,13,17,19]]))

figure3a

fig3b = Axis(figure3b[1,1], title="Cognitive score by quartile", xlabel="Age (years)")
leg3b = Legend(figure3b[1,2], [
    MarkerElement(color = colormap[1], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[16], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[3], marker = :circle, strokecolor = :black),
    ], ["lower 25%", "middle 50%", "upper 25%"], "Quartile", width = 150)
plot!(fig3b,
        collect(skipmissing(kids_metadata[hascog, :correctedAgeMonths])) ./ 12,
        collect(skipmissing(kids_metadata[hascog, :cogScore])), 
        color=categorical_colors(kids_metadata[hascog, :cogQuant], 1:3, colormap[[1,16,3]]))

figure3b

CairoMakie.save("figures/05_cogscore_age_test.svg", figure3a)
CairoMakie.save("figures/05_cogscore_age_quant.svg", figure3b)


#-
hasmullen = map(x-> x === "Mullen", kids_metadata.cogAssessment)
hasbayley = map(x-> x === "Bayleys", kids_metadata.cogAssessment)
haswppsi = map(x-> x === "WPPSI", kids_metadata.cogAssessment)

figure4 = Figure(resolution=(1200,800))
mullen = Axis(figure4[1,1], title = "age distribution for tests")

hist!(mullen, kids_metadata.correctedAgeDays[hasmullen], color=:blue)
hist!(mullen, kids_metadata.correctedAgeDays[hasbayley], color=:green)
hist!(mullen, kids_metadata.correctedAgeDays[haswppsi], color=:red)

extrema(kids_metadata.correctedAgeDays[hasmullen])
extrema(kids_metadata.correctedAgeDays[hasbayley])
extrema(kids_metadata.correctedAgeDays[haswppsi])
filter(row-> row.subject == 390, kids_metadata)

sum(hasbayley)
