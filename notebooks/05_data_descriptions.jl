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


figure1 = Figure(resolution=(900,900))

fig1a = Axis(figure1[1,1:2], title="Child ages", ylabel = "Count", xlabel="Age (Years)", xminorticksvisible=true, xminorgridvisible=true, xticks = 0:5:15, xminorticks=IntervalsBetween(5))
fig1b = Axis(figure1[2,1], title = "WHO BMI", ylabel = "BMI (z-score)", xlabel="Age (Years)", xminorticksvisible=true, xminorgridvisible=true, xticks = 0:1:5, xminorticks=IntervalsBetween(2))
fig1c = Axis(figure1[2,2], title = "Cognitive score", ylabel = "Score", xlabel="Age (Years)", xminorticksvisible=true, xminorgridvisible=true, xticks = 0:5:15, xminorticks=IntervalsBetween(5))
fig1c_leg = Legend(figure1[2,3], [
    MarkerElement(color = colormap[3], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[4], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[6], marker = :circle, strokecolor = :black),
    MarkerElement(color = colormap[9], marker = :circle, strokecolor = :black),
    ], ["Mullen", "Bayleys", "WPPSI", "WISC"], "Assessment")
# fig1i = Axis(figure1[5,1], title="Child ages 1-2", ylabel = "Count", xminorticksvisible=true, xminorgridvisible=true, xticks = 0:0.5:2, xminorticks=IntervalsBetween(6))

hist!(fig1a, kids_metadata.correctedAgeYears, color=:gray)

hasbmi = .!ismissing.(kids_metadata.WHO_zbmi)
scatter!(fig1b, collect(skipmissing(kids_metadata[hasbmi, :correctedAgeYears])), collect(skipmissing(kids_metadata[hasbmi, :WHO_zbmi])))

hascog = .!ismissing.(kids_metadata.cogScore)
scatter!(fig1c, collect(skipmissing(kids_metadata[hascog, :correctedAgeYears])), collect(skipmissing(kids_metadata[hascog, :cogScore])), 
        color=categorical_colors(kids_metadata[hascog, :cogAssessment], ["Mullen", "Bayleys", "WPPSI", "WISC"], colormap[[3,4,6,9]]))

tightlimits!(fig1a, Bottom())
figure1

ax1 = Axis(figure1[3,1], ylabel="First timepoint", xlabel="count")
let stk = stack(counts, [:single,:multiple])
    stk.order = map(x-> x== "single" ? 1 : 2, stk.variable)
    barplot!(ax1, stk.first_tp, stk.value, stack=stk.order, color=repeat(colormap[[2,8]], inner=4), direction=:x)
    ax1.yticks = (1:4, unique(stk.first_tp))
    leg = Legend(figure1[4,1], [
        MarkerElement(color = colormap[2], marker = :rect, strokecolor = :black),
        MarkerElement(color = colormap[8], marker = :rect, strokecolor = :black),
        ], ["Single", "Multiple"], "Number of timepoints",
        orientation=:horizontal, tellwidth=false, tellheight=true)
    tightlimits!(ax1, Left())
end 
ax2 = Axis(figure1[3,2], ylabel="First timepoint", xlabel="count")
let stk = stack(counts, [:infant, :early, :middle, :adolescent])
    stk.order = map(x-> x == "infant" ? 1 : 
                        x == "early"  ? 2 :
                        x == "middle" ? 3 :
                                        4, stk.variable)
    cs = ColorSchemes.seaborn_colorblind.colors
    barplot!(ax2, stk.first_tp, stk.value, stack=stk.order, color=repeat(cs[[1,2,3,5]], inner=4), direction=:x)
    hideydecorations!(ax2)
    leg = Legend(figure1[4,2], [
        MarkerElement(color = cs[1], marker = :rect, strokecolor = :black),
        MarkerElement(color = cs[2], marker = :rect, strokecolor = :black),
        MarkerElement(color = cs[3], marker = :rect, strokecolor = :black),
        MarkerElement(color = cs[5], marker = :rect, strokecolor = :black),
        ], ["infant", "early", "middle", "adolescent"], "Last timepoint",
        orientation=:horizontal, tellwidth=false, tellheight=true, nbanks=2)
    tightlimits!(ax2, Left())
end
        
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

floatshrinker(v,i,j) = v isa AbstractFloat ? round(v, digits=2) : v

kids_metadata.agePeriod = map(kids_metadata.correctedAgeDays) do a
    years = floor(Int, a / 365)
    years < 1 ? "Infancy" :
    years < 6 ? "Early childhood" :
    years < 12 ? "Middle childhood" :
                 "Adolecence"
end

@chain kids_metadata begin
    groupby(:agePeriod)
    combine(:correctedAgeYears => mean => "Mean", 
            :correctedAgeYears => median => "Median",
            :correctedAgeYears => minimum => "Minimum",
            :correctedAgeYears => maximum => "Maximum",
            nrow => "Total")
    sort(:agePeriod, lt=(x,y) -> x == "Infancy" ? true : 
                                 y == "Infancy" ? false :
                                 x == "Early childhood" ? true : 
                                 y == "Early childhood" ? false :
                                 x == "Middle childhood" ? true : 
                                 false
                                 )
    rename(:agePeriod=> " ")
    pretty_table(String, _; backend=:latex, nosubheader=true,
                label="tab:agestats", formatters=floatshrinker)
    print
end

@chain unique(kids_metadata, :subject) begin
    groupby(:simple_race)
    combine(nrow=> "Number")
    sort(:Number)
    rename(:simple_race=> :Race)
    pretty_table(String, _; backend=:latex, nosubheader=true, label="tab:race")
    print
end

@chain unique(has_bf, :subject) begin
    groupby(:breastfeeding)
    combine(nrow=> "Number")
    sort(:Number)
    rename(:breastfeeding=> "Liquid diet")
    pretty_table(String, _; backend=:latex, nosubheader=true, label="tab:breastfeeding")
    print
end

@chain unique(kids_metadata, :subject) begin
    groupby(:birthType)
    combine(nrow=> "Number")
    sort(:Number)
    rename(:birthType=> :Delivery)
    pretty_table(String, _; backend=:latex, nosubheader=true, label="tab:delivery")
    print
end
