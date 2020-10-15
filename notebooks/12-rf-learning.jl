include("../scripts/startup_loadpackages.jl")
include("../scripts/accessories.jl")
using MLJ


config = parsefile("Data.toml")
allmeta = CSV.File(config["tables"]["joined_metadata"], pool=false) |> DataFrame

filter!(:correctedAgeDays=> !ismissing, allmeta)
let samples = Set(sampleid.(uniquetimepoints(stoolsample.(allmeta.sample), takefirst=false, samplefilter=iskid)))
    filter!(:sample => s-> s ∈ samples, allmeta)  
end

allmeta.ageLabel = map(eachrow(allmeta)) do row
    ismissing(row.correctedAgeDays) && error("No age for $(row.sample)")
    row.correctedAgeDays < 365 && return "1 and under"
    row.correctedAgeDays < 365*2 && return "1 to 2"
    return "2 and over"
end

filter!([:cogScore, :braintotal]=> (cs, bt)-> any(!ismissing, (cs, bt)), allmeta)

species = widen2comm(taxonomic_profiles(filefilter=f-> sampleid(stoolsample(basename(f))) in allmeta.sample)...)
filter!(:sample => s-> s ∈ samplenames(species), allmeta)
filter!(:cogScore => !ismissing, allmeta)

species = view(species, sites=allmeta.sample)
@assert allmeta.sample == samplenames(species)


(l, m, u) = quantile(skipmissing(allmeta.cogScore), [0.25,0.5,0.75])
DataFrames.transform!(allmeta, :cogScore => ByRow(s-> ismissing(s) ? missing :
                            s < l ? "bottom" : 
                            s < m ? "lower-mid" : 
                            s < u ? "upper-mid" :
                        "upper") => :cogQuartile).cogQuartile

allmeta.cogQuartile = categorical(allmeta.cogQuartile)
levels!(allmeta.cogQuartile, ["bottom", "lower-mid", "upper-mid", "upper"])

sigbugs = CSV.File(joinpath(config["output"]["tables"], "oldkidsquartiletests.csv")) |> DataFrame
sigbugs = Symbol.(sigbugs[sigbugs.qvalue .< 0.2, :species])

for sp in String.(sigbugs)
    v = vec(occurrences(view(species, sites=allmeta.sample, species=[sp])))
    allmeta[:, sp] = v
end


oldkids = view(allmeta, allmeta.correctedAgeDays .> 365, :)

ukids = let samples = Set(sampleid.(uniquetimepoints(oldkids.sample, takefirst=true, samplefilter=iskid)))
    findall(s-> in(s, samples), oldkids.sample)
end


y, X = unpack(oldkids, ==(:cogQuartile), col-> col in sigbugs)

tree_model = MLJ.@load DecisionTreeClassifier verbosity=1

ev = MLJ.evaluate(tree_model, X, y,
                resampling=CV(shuffle=true), measure=cross_entropy, verbosity=0)


tree = machine(tree_model, X, y)

MLJ.fit!(tree, rows=ukids)
yhat = MLJ.predict(tree, X[Not(ukids),:])
cross_entropy(yhat, y[Not(ukids)]) |> mean
MLJ.confusion_matrix(mode.(yhat), y[Not(ukids)])

