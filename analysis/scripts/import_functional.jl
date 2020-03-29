using SQLite
using CSV
using Pkg.TOML: parsefile
using ECHOAnalysis

allsamples = readdir(config["files"]["rawfastq"]["path"])
filter!(isstoolsample, allsamples)
allsamples = unique(stoolsample.(allsamples))
samples = uniquetimepoints(allsamples, takefirst=false)

funcdb = SQLite.DB(config["sqlite"]["uniref90"]["path"])
add_functional_profiles(funcdb, config["files"]["biobakery"]["path"], kind="genefamilies_relab", replace=true, samples=samples)
kodb = SQLite.DB(config["sqlite"]["ko"]["path"])
add_functional_profiles(kodb, config["files"]["biobakery"]["path"], kind="ko_names_relab", replace=true, samples=samples)
pfamdb = SQLite.DB(config["sqlite"]["pfam"]["path"])
add_functional_profiles(pfamdb, config["files"]["biobakery"]["path"], kind="pfam_names_relab", replace=true, samples=samples)
ecdb = SQLite.DB(config["sqlite"]["ec"]["path"])
add_functional_profiles(ecdb, config["files"]["biobakery"]["path"], kind="ec_names_relab", replace=true, samples=samples)
