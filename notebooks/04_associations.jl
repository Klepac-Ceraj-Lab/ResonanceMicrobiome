using ResonanceMicrobiome
using BiobakeryUtils
using CairoMakie
using AbstractPlotting.ColorSchemes

#-

all_species = taxonomic_profiles(:species)  
all_unirefs = functional_profiles(:uniref)
all_metadata = resonance_metadata(name.(samples(all_species)))