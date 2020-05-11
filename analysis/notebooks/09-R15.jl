include("../scripts/startup_loadpackages.jl")
@load "analysis/figures/assets/metadata.jld2" allmeta ubothmeta ukidsmeta allkidsmeta allmoms allkids umoms ukids oldkids uboth
@load "analysis/figures/assets/taxa.jld2" species speciesmds kidsspeciesmds kidsspeciesmdsaxes
@load "analysis/figures/assets/unirefs.jld2" unirefaccessorymds unirefaccessorymdsaxes kidsunirefaccessorymds kidsunirefaccessorymdsaxes
@load "analysis/figures/assets/otherfunctions.jld2" kos kosdiffs kosdm ecs ecsdm pfams pfamsdiffs pfamsdm
@load "analysis/figures/assets/permanovas.jld2" r2 r2m qa allpermanovas species_permanovas unirefaccessory_permanovas kos_permanovas pfams_permanovas
@load "analysis/figures/assets/fsea.jld2" allfsea mdcors
@load "analysis/figures/assets/difs.jld2" speciesdiffs unirefaccessorydiffs kosdiffs pfamsdiffs

ebf = findall(x-> !ismissing(x) && x == "exclussive breast", allkidsmeta.breastfeeding)
eff = findall(x-> !ismissing(x) && x == "exclussive formula", allkidsmeta.breastfeeding)

allkidsspecies = view(species, sites=allkids)

feature = view(allkidsspecies, species=[1]) |> occurrences
feature_bf = feature[ebf]
feature_ff = feature[eff]
SignedRankTest(feature_bf, feature_ffs)

test = CSV.File(config["tables"]["joined_metadata"], pool=false) |> DataFrame
