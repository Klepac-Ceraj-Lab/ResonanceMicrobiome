using CSV
using DataFrames
using AbstractPlotting
using AbstractPlotting.MakieLayout
using CairoMakie
using StatsMakie

include(joinpath("..", "scripts", "startup_loadpackages.jl"))
include(joinpath("..", "scripts", "accessories.jl"))
config = parsefile("Data.toml")
allmeta = CSV.File(config["tables"]["joined_metadata"], pool=false) |> DataFrame

filter!(:correctedAgeDays=> !ismissing, allmeta)

hires = CSV.File(joinpath("data", "brain", "microbiome_corticalVolumes_bambam.csv")) |> DataFrame
hires2 = CSV.File(joinpath("data", "brain", "microbiome_corticalVolumes_bambam_sean.csv")) |> DataFrame
rename!(hires, Dict(:ID=>:subject, :Timepoint=>:timepoint))
rename!(hires2, Dict(:ID=>:subject, :Timepoint=>:timepoint))


hr2_samples = resolve_letter_timepoint.(string.(hires2.subject))
hires2.subject = subject.(hr2_samples)
hires2.timepoint = timepoint.(hr2_samples)

# don't want to replicate :age column
select!(hires, Not(:Age))
@assert names(hires) == names(hires2)
hires = vcat(hires, hires2)
unique!(hires, [:subject,:timepoint])

# There are a lot of individual brain regions that are separated in this table,
# and the right and left hemispheres are distinguished.
# For the most part, we're not going to need this level of specificity,
# but we can group individual brain regions
# and combine left / right hemispheres.
# I'll also make a column with the total brain volume for later normalization.

mapping = CSV.File("data/brain/brain_region_key.csv") |>DataFrame
hires.hires_total = [sum(row[3:end]) for row in eachrow(hires)]

cols_seen = let ns = lowercase.(String.(names(hires)))
    cols_seen = Int[]
    by(mapping, :region) do region
        fs = lowercase.(region.feature)
        cols = findall(n-> any(f-> occursin(f, n), fs), ns)
        append!(cols_seen, cols)
        hires[!, Symbol(first(region.region))] = [sum(row[cols]) for row in eachrow(hires)]
        true # need something for the `by`
    end
    cols_seen
end

select!(hires, ["subject", "timepoint", "cerebellar", "neocortical", "subcortical", "limbic", "hires_total"])
rename!(hires, ["cerebellar"=> "cerebellum_old", "neocortical"=> "neocortex_old", "limbic"=> "limbic_old", "subcortical"=> "subcortex_old"])

allmeta = leftjoin(allmeta, hires, on=["subject", "timepoint"])
nomissing = filter(allmeta) do row
    !any(ismissing, [
        row."cerebellum",
        row."cerebellum_old",
        row."limbic",
        row."limbic_old",
        row."neocortex",
        row."neocortex_old",
        row."subcortex",
        row."cerebellum_old",
        row."correctedAgeDays"
    ])
end

disallowmissing!(nomissing, [ "neocortex", "neocortex_old", "subcortex", "subcortex_old", "cerebellum", "cerebellum_old", "limbic", "limbic_old", "braintotal", "hires_total"])

## 
res = (1000, 1500)
scene, layout = layoutscene(resolution=res)

ax1 = layout[1,1] = LAxis(scene, title="neocortex", xlabel="new", ylabel="old")
ax2 = layout[1,2] = LAxis(scene, title="subcortex", xlabel="new", ylabel="old")
ax3 = layout[2,1] = LAxis(scene, title="cerebellum", xlabel="new", ylabel="old")
ax4 = layout[2,2] = LAxis(scene, title="limbic", xlabel="new", ylabel="old")
ax5 = layout[3,1] = LAxis(scene, title="totals", xlabel="new (icv)", ylabel="old (sum)")

function plotcor!(ax, x1, x2)
    (low, high) = extrema([x1; x2])
    plot!(ax, [low, high], [low, high])
    xlims!(ax, low, high),
    ylims!(ax, low, high)
end

scatter!(ax1, nomissing.neocortex, nomissing.neocortex_old ./ 1000)
plotcor!(ax1, nomissing.neocortex, nomissing.neocortex_old ./ 1000)
scatter!(ax2, nomissing.subcortex, nomissing.subcortex_old)
plotcor!(ax2, nomissing.subcortex, nomissing.subcortex_old)
scatter!(ax3, nomissing.cerebellum, nomissing.cerebellum_old)
plotcor!(ax3, nomissing.cerebellum, nomissing.cerebellum_old)
scatter!(ax4, nomissing.limbic, nomissing.limbic_old)
plotcor!(ax4, nomissing.limbic, nomissing.limbic_old)
scatter!(ax5, nomissing.braintotal, nomissing.hires_total)
plotcor!(ax5, nomissing.braintotal, nomissing.hires_total)

scene
save("/home/kevin/Desktop/brain_descrepancies1_fixneo.pdf", scene, resolution=res);

##
