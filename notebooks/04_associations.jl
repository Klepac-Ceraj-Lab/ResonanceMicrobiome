using ResonanceMicrobiome
using BiobakeryUtils
using CairoMakie
using AbstractPlotting.ColorSchemes



#-

all_species = taxonomic_profiles(:species)
all_species = all_species[:, map(s->!occursin(r"_\d+E_", s), samplenames(all_species))]
all_unirefs = functional_profiles(:uniref)
all_metadata = resonance_metadata(name.(samples(all_species)))
unique!(all_metadata, [:subject, :timepoint])
has_age = all_metadata[.!ismissing.(all_metadata.correctedAgeDays), :]
#-

fig = Figure();

fa = fig[1,1] = Axis(fig, xlabel="Age (Years)", ylabel="Count", title="All")
hist!(fa,  has_age.correctedAgeDays ./ 365, color=:grey)

fb = fig[1,2] = Axis(fig, xlabel="Age (Years)", ylabel="Count", title="Under 2")
hist!(fb,  filter(row-> row.correctedAgeDays < 365*2, has_age).correctedAgeDays ./ 365, color=:grey)

fig

dpm = 365 / 12

for i in 1:24
    inter = filter(row-> dpm * (i-1) < row.correctedAgeDays < dpm * i, has_age)
    count = nrow(inter)
    println("Between $(i-1) and $i months: $count")
end

threeseven = filter(row-> dpm*3 < row.correctedAgeDays < dpm*7, has_age).subject
foursix = filter(row-> dpm*4 < row.correctedAgeDays < dpm*6, has_age).subject
fiveseven = filter(row-> dpm*5 < row.correctedAgeDays < dpm*7, has_age).subject

unique(filter(row-> 
    (row.subject in fiveseven) &&
    11*dpm < row.correctedAgeDays < 13*dpm, has_age), :subject) |> nrow
unique(filter(row-> 
    (row.subject in fiveseven) &&
    11*dpm < row.correctedAgeDays < =*dpm, has_age), :subject) |> nrow
unique(filter(row-> 
    (row.subject in fiveseven) &&
    11*dpm < row.correctedAgeDays < 13*dpm, has_age), :subject) |> nrow


unique(filter(row-> row.subject in fiveseven, has_age), :subject)
