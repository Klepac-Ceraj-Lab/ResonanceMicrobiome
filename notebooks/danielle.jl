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
using SparseArrays

colormap = ColorSchemes.tab20.colors
all_metadata = CSV.read("../danielle-thesis/clinical_metadata.csv", DataFrame)
rename!(all_metadata, :subjectID=>:subject)

filter!(:correctedAgeDays=> !ismissing, all_metadata)
println.(names(all_metadata))
#-
mgx_genera = CSV.read("/home/kevin/repos/danielle-thesis/paper-taxonomic-levels/resonance_mgx_genus.csv", DataFrame)
mgx_genera = CommunityProfile(sparse(Matrix(mgx_genera[!, 2:end])), Taxon.(replace.(mgx_genera.taxname, Ref("g__"=>"")),:genus), MicrobiomeSample.(replace.(names(mgx_genera)[2:end], Ref("-"=>"_"))))
amp_genera = CSV.read("/home/kevin/repos/danielle-thesis/paper-taxonomic-levels/resonance_dada2_genera.csv", DataFrame)
amp_genera = CommunityProfile(sparse(Matrix(amp_genera[!, 2:end])), Taxon.(replace.(amp_genera.genus, Ref("g__"=>"")),:genus), MicrobiomeSample.(replace.(names(amp_genera)[2:end], Ref("-"=>"_"))))

shared_samples = filter(s-> !startswith(s, "M"), name.(samples(mgx_genera)) ∩ name.(samples(amp_genera)))
sharedmeta = DataFrame(sample = shared_samples)
sharedmeta.subject = map(s-> parse(Int, match(r"C(\d+)", s).captures[1]), shared_samples)
sharedmeta.timepoint = map(s-> parse(Int, match(r"_(\d+)[FE]_", s).captures[1]), shared_samples)
sharedmeta = leftjoin(sharedmeta, select(all_metadata, [:subject,:timepoint, :correctedAgeDays]), on = [:subject, :timepoint])
filter!(:correctedAgeDays=> !ismissing, sharedmeta)

mgx_genera = mgx_genera[:, sharedmeta.sample]
amp_genera = amp_genera[:, sharedmeta.sample]

sharedmeta.ageMonths = sharedmeta.correctedAgeDays ./ 365

mgx_15_genera = mgx_genera[:, filter(:ageMonths => a-> a <= 15, sharedmeta).sample]    
mgx_30_genera = mgx_genera[:, filter(:ageMonths => a-> 15 < a <= 30, sharedmeta).sample]    
mgx_o30_genera = mgx_genera[:, filter(:ageMonths => a-> a > 30, sharedmeta).sample]    
amp_15_genera = amp_genera[:, filter(:ageMonths => a-> a <= 15, sharedmeta).sample]    
amp_30_genera = amp_genera[:, filter(:ageMonths => a-> 15 < a <= 30, sharedmeta).sample]    
amp_o30_genera = amp_genera[:, filter(:ageMonths => a-> a > 30, sharedmeta).sample]    


mgx_15_bc = braycurtis(mgx_15_genera)
mgx_30_bc = braycurtis(mgx_30_genera)
mgx_o30_bc = braycurtis(mgx_o30_genera)
amp_15_bc = braycurtis(amp_15_genera)
amp_30_bc = braycurtis(amp_30_genera)
amp_o30_bc = braycurtis(amp_o30_genera)

mgx_15_clust = hclust(mgx_15_bc, linkage=:complete, branchorder=:optimal)
mgx_30_clust = hclust(mgx_30_bc, linkage=:complete, branchorder=:optimal)
mgx_o30_clust = hclust(mgx_o30_bc, linkage=:complete, branchorder=:optimal)
amp_15_clust = hclust(amp_15_bc, linkage=:complete, branchorder=:optimal)
amp_30_clust = hclust(amp_30_bc, linkage=:complete, branchorder=:optimal)
amp_o30_clust = hclust(amp_o30_bc, linkage=:complete, branchorder=:optimal)


function topx(cp, n=10)
    totals = vec(featuretotals(cp))
    rows = partialsortperm(totals, 1:n, rev=true)
    top = cp[rows, :]
    other = sum(abundances(cp[collect(1:nfeatures(cp))[Not(rows)], :]), dims=1)

    return (names = vcat(featurenames(top), ["other"]), abundances = vcat(abundances(top), other))
end

allnames = setdiff(union([topx(cp).names for cp in (mgx_15_genera, 
                                                    mgx_30_genera, 
                                                    mgx_o30_genera, 
                                                    amp_15_genera, 
                                                    amp_30_genera, 
                                                    amp_o30_genera)]...), Set(["other"]))
name_dict = Dict(n => colormap[i] for (i, n) in enumerate(allnames))
name_dict["other"] = ColorSchemes.Greys_3.colors[1]

function plottopn!(fig, layout, axrow, cp, clust, n, title, markerspecies = String[]; kwargs...)
    topnames, topabund = topx(cp, n)
    topabund = topabund[:, clust.order]
    ax = layout[axrow,1] = Axis(fig; title)
    append!(markerspecies, [topnames[i] for i in 1:n+1])
    
    for i in (n+1):-1:1
        v = vec(sum(topabund[1:i, :], dims=1))
        barplot!(ax, 1:nsamples(cp), v, color=name_dict[topnames[i]])
    end

    tightlimits!(ax)
    hidexdecorations!(ax)
    fig, markerspecies
end

sublayout = GridLayout()
figure1[1:2, 4:6] = sublayout
(_, markerspecies) = plottopn!(figure1, sublayout, 1, u1_spec, u1clust, 10, "Top 10 species, kids under 1 yo")
plottopn!(figure1, sublayout, 2, mid_spec, midclust, 10, "Top 10 species, kids 1-2 yo", markerspecies)
plottopn!(figure1, sublayout, 3, o2_spec, o2clust, 10, "Top 10 species, kids over 2 yo", markerspecies)
unique!(markerspecies)
sort!(markerspecies)
fig1c_leg = sublayout[4, 1] = Legend(figure1, [MarkerElement(color = name_dict[m], marker = :rect, strokecolor = :black) for m in markerspecies], [m for m in markerspecies], tellheight=true, nbanks = 2)

CairoMakie.save("figures/03_taxonomic_profiles.pdf", figure1)
figure1
