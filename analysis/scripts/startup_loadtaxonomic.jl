### Load taxonomic profiles from SQLite database
taxdb = SQLite.DB(config["sqlite"]["taxa"]["path"])
### Function in ECHOAnalysis package
species = sqlprofile(taxdb, tablename="taxa", kind="species")
### Keep only samples in the metadata
species = view(species, sites=allmeta.sample) |> copy
### Total sum scaling - function in Microbiome
relativeabundance!(species)
