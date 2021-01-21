module ResonanceMicrobiome

export
    @datadep_str,
    resonance_metadata

using Reexport 
@reexport using Microbiome
@reexport using DataFrames
@reexport using CSV

using Arrow
using DataDeps
import DataDeps: @datadep_str
using LoggingExtras
using CairoMakie
using Airtable


include("data.jl")


end # module