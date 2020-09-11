using ECHOAnalysis
using DataFrames
using CSV
import Pkg.TOML: parsefile 

config = parsefile("Data.toml")
allmeta = CSV.File(config["tables"]["joined_metadata"], pool=false) |> DataFrame
allmeta.batch

taxa = readlines(joinpath(config["output"]["other"], "echo_taxa.txt"))

