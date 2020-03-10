using SQLite
using CSV
using Pkg.TOML: parsefile
using ECHOAnalysis

config = parsefile("data/data.toml")

taxdb = SQLite.DB(config["sqlite"]["taxa"]["path"])
add_taxonomic_profiles(taxdb, config["files"]["biobakery"]["path"], replace=true)
