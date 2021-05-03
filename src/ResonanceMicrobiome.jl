module ResonanceMicrobiome

export
    @datadep_str,
    resonance_metadata,
    taxonomic_profiles,
    functional_profiles,
    sample_filter,
    sample_filter!,
    clinical_data_descriptions

export categorical_colors,
       mds_format,
       mds_percent

using Reexport 
@reexport using Microbiome
@reexport using DataFrames
@reexport using CSV

using SparseArrays
using Arrow
using DataDeps
import DataDeps: @datadep_str
using LoggingExtras
using CairoMakie
using Airtable
using BiobakeryUtils

include("data_sources.jl")
include("data_import.jl")
include("plot_styling.jl")

end # module