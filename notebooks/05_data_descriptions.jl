using ResonanceMicrobiome
using BiobakeryUtils
using CairoMakie
using StatsBase
using SplitApplyPlot
using SplitApplyPlot: CategoricalScale
using AbstractPlotting.ColorSchemes
using Statistics
using CategoricalArrays

CairoMakie.activate!(type="svg")
colormap = ColorSchemes.tab20.colors

SplitApplyPlot.default_axis(::Type{Axis}) = NamedTuple()
SplitApplyPlot.default_axis(::Type{Axis3}) = NamedTuple()
SplitApplyPlot.default_styles() = NamedTuple()

#-

all_metadata = resonance_metadata()
all_metadata.correctedAgeYears = all_metadata.correctedAgeDays ./ 365
# all_metadata.breastfeeding = categorical(all_metadata.breastfeeding)
# levels!(all_metadata.breastfeeding, ["exclusive formula", "exclusive breast", "mixed"])

kids_metadata = filter(row-> !ismissing(row.correctedAgeDays), all_metadata)
has_bf = filter(row-> !ismissing(row.breastfeeding), all_metadata)

unique(kids_metadata.subject)

#- longitudinal thing

subj = groupby(kids_metadata, :subject)
counts = DataFrame(first_tp = ["infant", "early", "middle", "adolescent"],
                    single = zeros(Int, 4),
                    multiple = zeros(Int, 4),
                    infant = zeros(Int, 4),
                    early = zeros(Int, 4),
                    middle = zeros(Int, 4),
                    adolescent = zeros(Int, 4))
for group in subj
    subdf = sort(group, :correctedAgeDays)
    started = ceil(Int, subdf[1,:correctedAgeDays] / 365)
    row = started <= 1  ? 1 : # infant
          started <= 6  ? 2 : # early childhood
          started <= 12 ? 3 : # late childhood
                          4   # adolescent
    nrow(subdf) > 1 ? counts[row, 3] += 1 : counts[row, 2] += 1
    for a in subdf.correctedAgeDays
        years = ceil(Int, a / 365)
        col = years <= 1  ? 4 : # infant
              years <= 6  ? 5 : # early childhood
              years <= 12 ? 6 : # late childhood
                            7   # adolescent
        counts[row, col] += 1
    end
end


figure1 = Figure(resolution=(1200,1600))

fig1a = Axis(figure1[1,1], title="Child ages", ylabel = "Count", xlabel="Age (Years)", xminorticksvisible=true, xminorgridvisible=true, xticks = 0:5:15, xminorticks=IntervalsBetween(5))
fig1b = Axis(figure1[1,2], title="Child ages (0-2)", ylabel = "Count", xlabel="Age (months)", xminorticksvisible=true, xminorgridvisible=true, xticks = 0:6:24, xminorticks=IntervalsBetween(6))
fig1c = Axis(figure1[2,1], ylabel = "Count")
fig1c_leg = Legend(figure1[2,3], [
    MarkerElement(color = colormap[1], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[7], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[9], marker = :circle, strokecolor = :black),
    ], ["Breastmilk", "Formula", "Mixed"], "Diet")

fig1e = Axis(figure1[3,1], title = "Race / ethnicity", ylabel = "Count")
fig1f = Axis(figure1[3,2], title = "SES", ylabel = "Count")
fig1g = Axis(figure1[4,1], title = "WHO BMI", ylabel = "BMI (z-score)", xlabel="Age (Years)", xminorticksvisible=true, xminorgridvisible=true, xticks = 0:1:5, xminorticks=IntervalsBetween(2))
fig1h = Axis(figure1[4,2], title = "Cognitive score", ylabel = "Score", xlabel="Age (Years)", xminorticksvisible=true, xminorgridvisible=true, xticks = 0:5:15, xminorticks=IntervalsBetween(5))
fig1h_leg = Legend(figure1[4,3], [
    MarkerElement(color = colormap[10], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[13], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[17], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[19], marker = :circle, strokecolor = :black),
    ], ["Mullen", "Bayleys", "WPPSI", "WISC"], "Assessment")
# fig1i = Axis(figure1[5,1], title="Child ages 1-2", ylabel = "Count", xminorticksvisible=true, xminorgridvisible=true, xticks = 0:0.5:2, xminorticks=IntervalsBetween(6))

fig1j = Axis(figure1[5,2], ylabel = "Count")
fig1j_leg = Legend(figure1[5,3], [
    MarkerElement(color = colormap[end-2], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[end], marker = :circle, strokecolor = :black),
    ], ["Cesarean", "Vaginal"], "Delivery method")
    


hist!(fig1a, kids_metadata.correctedAgeYears, color=:gray)
hist!(fig1b, filter(:correctedAgeYears=> <(2), kids_metadata).correctedAgeYears .* 12, color=:gray)

barplot!(fig1c, [1,2,3], [
    count(x-> x === ("exclusive breast"), all_metadata.breastfeeding),
                count(x-> x === ("exclusive formula"), all_metadata.breastfeeding),
                count(x-> x === ("mixed"), all_metadata.breastfeeding)],
                color = colormap[[1,7,9]])
fig1c.xticks = ([1,2,3], ["Formula", "Breastmilk", "Mixed"])
fig1c.xticklabelsize = 20

let specs = data(has_bf) * 
        mapping(:correctedAgeYears=> "Age (Years)", 
            stack=:breastfeeding, 
            color=:breastfeeding) *
        histogram(bins=20)
    test = draw!(figure1[2,2], specs; palettes=(;color=colormap[[1,7,9]]))
    tightlimits!(test[1].axis)
end

let fig = figure1[5,1] = GridLayout()
    ax1 = fig[1,1] = Axis(figure1, ylabel="First timepoint", xlabel="count")
    stk = stack(counts, [:single,:multiple])
    stk.order = map(x-> x== "single" ? 1 : 2, stk.variable)
    barplot!(ax1, stk.first_tp, stk.value, stack=stk.order, color=repeat(colormap[[9,14]], inner=4), direction=:x)
    ax1.yticks = (1:4, unique(stk.first_tp))
    leg = fig[2,1] = Legend(figure1, [
        MarkerElement(color = colormap[9], marker = :rect, strokecolor = :black),
        MarkerElement(color = colormap[14], marker = :rect, strokecolor = :black),
        ], ["Single", "Multiple"], "Number of timepoints",
        orientation=:horizontal, tellwidth=false, tellheight=true)
        tightlimits!(ax1, Left())
 
    ax2 = fig[1,2] = Axis(figure1, ylabel="First timepoint", xlabel="count")
    stk = stack(counts, [:infant, :early, :middle, :adolescent])
    stk.order = map(x-> x == "infant" ? 1 : 
                        x == "early"  ? 2 :
                        x == "middle" ? 3 :
                                        4, stk.variable)
    barplot!(ax2, stk.first_tp, stk.value, stack=stk.order, color=repeat(colormap[[11,12,15,16]], inner=4), direction=:x)
    hideydecorations!(ax2)
    leg = fig[2,2] = Legend(figure1, [
        MarkerElement(color = colormap[11], marker = :rect, strokecolor = :black),
        MarkerElement(color = colormap[12], marker = :rect, strokecolor = :black),
        MarkerElement(color = colormap[15], marker = :rect, strokecolor = :black),
        MarkerElement(color = colormap[16], marker = :rect, strokecolor = :black),
        ], ["infant", "early", "middle", "adolescent"], "Last timepoint",
        orientation=:horizontal, tellwidth=false, tellheight=true, nbanks=2)
    tightlimits!(ax2, Left())
end

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
scatter!(fig1g, collect(skipmissing(kids_metadata[hasbmi, :correctedAgeYears])), collect(skipmissing(kids_metadata[hasbmi, :WHO_zbmi])))

hascog = .!ismissing.(kids_metadata.cogScore)
scatter!(fig1h, collect(skipmissing(kids_metadata[hascog, :correctedAgeYears])), collect(skipmissing(kids_metadata[hascog, :cogScore])), 
        color=categorical_colors(kids_metadata[hascog, :cogAssessment], ["Mullen", "Bayleys", "WPPSI", "WISC"], colormap[[10,13,17,19]]))


barplot!(fig1j, [1,2], [
    count(x-> x === ("Cesarean"), unique(all_metadata, :subject).birthType),
    count(x-> x === ("Vaginal"), unique(all_metadata, :subject).birthType)],
    color = colormap[[end-2,end]])
fig1j.xticks = ([1,2], ["Cesarean", "Vaginal"])
fig1j.xticklabelsize = 20

tightlimits!.([fig1a, fig1b, fig1c, fig1e, fig1f, fig1j])
figure1
CairoMakie.save("figures/05_data_summaries.pdf", figure1)

#- 

has_16S = filter(row-> !ismissing(row."16S_batch"), kids_metadata)
figure2 = Figure(resolution=(1200,800))

age16S = Axis(figure2[1,1], title = "age distribution for 16S")
hist!(age16S, has_16S.correctedAgeYears)
figure2

describe(has_16S.correctedAgeYears)
count(a -> 4 < a < 6, has_16S.correctedAgeYears)
count(a -> 11 < a < 13, has_16S.correctedAgeYears)
filter(row -> 4 < row.correctedAgeYears < 6, has_16S).subject ∩ filter(row -> 11 < row.correctedAgeYears < 13, has_16S).subject

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
        collect(skipmissing(kids_metadata[hascog, :correctedAgeYears])) ./ 12,
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
        collect(skipmissing(kids_metadata[hascog, :correctedAgeYears])) ./ 12,
        collect(skipmissing(kids_metadata[hascog, :cogScore])), 
        color=categorical_colors(kids_metadata[hascog, :cogQuant], 1:3, colormap[[1,16,3]]))

figure3b

CairoMakie.save("figures/05_cogscore_age_test.pdf", figure3a)
CairoMakie.save("figures/05_cogscore_age_quant.pdf", figure3b)

#-
hasmullen = map(x-> x === "Mullen", kids_metadata.cogAssessment)
hasbayley = map(x-> x === "Bayleys", kids_metadata.cogAssessment)
haswisc = map(x-> x === "WISC", kids_metadata.cogAssessment)
haswppsi = map(x-> x === "WPPSI", kids_metadata.cogAssessment)

figure4 = Figure(resolution=(1200,800))
mullen = Axis(figure4[1,1], title = "age distribution for tests")

hist!(mullen, kids_metadata.correctedAgeDays[hasmullen], color=:blue)
hist!(mullen, kids_metadata.correctedAgeDays[hasbayley], color=:green)
hist!(mullen, kids_metadata.correctedAgeDays[haswppsi], color=:red)
hist!(mullen, kids_metadata.correctedAgeDays[haswisc], color=:purple)
figure4
extrema(kids_metadata.correctedAgeDays[hasmullen])
extrema(kids_metadata.correctedAgeDays[hasbayley])
extrema(kids_metadata.correctedAgeDays[haswppsi])
filter(row-> row.subject == 390, kids_metadata)

sum(hasbayley)

#- Numbers


using Chain
using Statistics
using DataFrames.PrettyTables

@chain kids_metadata begin
    groupby(:ageLabel)
    combine(:correctedAgeYears => mean => "Mean", 
            :correctedAgeYears => median => "Median",
            :correctedAgeYears => minimum => "Minimum",
            :correctedAgeYears => maximum => "Maximum",
            nrow => "Total")
    sort(:ageLabel)
    rename(:ageLabel=> " ")
    pretty_table(String, _; backend=:latex, nosubheader=true)
    print
end

@chain kids_metadata begin
    groupby(:simple_race)
    combine(nrow=> "Number")
    sort(:Number)
    rename(:simple_race=> :Race)
    pretty_table(String, _; backend=:latex, nosubheader=true)
    print
end

@chain unique(kids_metadata) begin
    groupby(:breastfeeding)
    combine(nrow=> "Number")
    sort(:Number)
    rename(:breastfeeding=> "Liquid diet")
    pretty_table(String, _; backend=:latex, nosubheader=true)
    print
end