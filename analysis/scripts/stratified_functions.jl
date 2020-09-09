import StatsMakie.StructArrays: uniquesorted

include("accessories.jl")
nko = get_neuroactive_kos()
# replace this with `select!()` once it's added https://github.com/JuliaData/DataFrames.jl/pull/2080

function longfromcomm(cm; stratified=false)
    features = featurenames(cm)
    samples  = samplenames(cm)
    if stratified
        taxa = Union{Nothing,String}[]
        for (i, f) in enumerate(features)
            m = match(r"^(.+?)\|?(?:(?:g__([\w]+)\.s__([\w]+))?|(unclassified))$", f)
            (ko, genus, species, unclass) = m.captures
            if isnothing(genus) && isnothing(unclass)
                ko == f || error("Weird function $f")
                taxon = nothing
            else
                isnothing(unclass) ? taxon = "$species" : taxon = "unclassified"
            end
            ko = split(ko, ':')[1]
            features[i] = ko
            push!(taxa, taxon)
        end
    else
        taxa = fill("unlassified", length(features))
    end
    occ = occurrences(cm)

    df = DataFrame((sample=samples[j], func=features[i], taxon=taxa[i], abundance=occ[i,j])
                        for i in eachindex(features)
                        for j in eachindex(samples) if occ[i,j] > 0)
    return df
end

koslong = longfromcomm(stratkos, stratified=true)
rename!(koslong, :func=>:ko)
filter!(row->row.abundance > 0, koslong)

# koslong = hcat(koslong, DataFrame(map(eachrow(koslong)) do row
#     m = match(r"^(.+?)\|?(?:(?:g__([\w]+)\.s__([\w]+))?|(unclassified))$", row.function)
#     isnothing(m) && @show row.function
#     (ko, genus, species, unclass) = m.captures
#     if isnothing(genus) && isnothing(unclass)
#         ko == row.function || error("Weird function $(row.function)")
#         taxon = nothing
#     else
#         isnothing(unclass) ? taxon = "$species" : taxon = "unclassified"
#     end
#     ko = split(ko, ':')[1]
#     (ko=ko, taxon=taxon)
#     end
# ))

filter!(row-> !ismissing(row.ageLabel) && row.ageLabel != "mom", allmeta)
filter!(row-> !ismissing(row.cogScore), allmeta)

let samples = Set(allmeta.sample)
    filter!(row-> row.sample in samples, koslong)
end

agelabels = dictionary(sample=>age for (sample, age) in eachrow(allmeta[!, [:sample, :ageLabel]]))
koslong.agelabel = [agelabels[s] for s in koslong.sample]
cogscores = dictionary(sample=>score for (sample, score) in eachrow(allmeta[!, [:sample, :cogScore]]))
koslong.cogScore = [cogscores[s] for s in koslong.sample]

##

gluts = let geneset = Set(nko["Glutamate synthesis"])
    filt = filter(row-> row.ko in geneset, koslong)
    samples = by(filt, [:sample, :taxon]) do sample
        (total=sum(sample.abundance), agelabel=first(sample).agelabel, cogScore=first(sample).cogScore)
    end

    sort!(samples, [:agelabel, :sample])
    totals = filter(row-> isnothing(row.taxon), samples)
    filter!(row-> !isnothing(row.taxon), samples)
    totals.ord = invperm(sortperm(totals, [:agelabel, :total]))

    spectotals = by(samples, :taxon, spectotal = :total => sum)
    sort!(spectotals, :spectotal)
    topspec = Set(last(spectotals, 5).taxon)

    by(samples, :sample) do df
        df = filter(row-> !in(row.taxon, topspec), df)
        if size(df, 1) == 0
            1
        else
            other = sum(df.total)
            sample = first(df).sample
            push!(samples, (sample=sample, taxon="other", total=other, agelabel=agelabels[sample], cogScore=cogscores[sample]))
            1
        end
    end
    push!(topspec, "other")
    filter!(row-> row.taxon in topspec, samples)


    for s in unique(samples.sample)
        missing_tax = setdiff(topspec, samples[samples.sample .== s, :taxon])
        for t in missing_tax
            push!(samples, (sample=s, taxon=t, total=0., agelabel=agelabels[s], cogScore=cogscores[s]))
        end
    end

    xs = dictionary(s=> i for (s,i) in eachrow(totals[!,[:sample,:ord]]))
    samples.x = [xs[s] for s in samples.sample]
    sort!(samples, :x)
    samples
end

glutd = let geneset = Set(nko["Glutamate degradation"])
    filt = filter(row-> row.ko in geneset, koslong)
    samples = by(filt, [:sample, :taxon]) do sample
        (total=sum(sample.abundance), agelabel=first(sample).agelabel)
    end
    @info names(samples)
    sort!(samples, [:agelabel, :sample])
    totals = filter(row-> isnothing(row.taxon), samples)
    filter!(row-> !isnothing(row.taxon), samples)
    totals.ord = invperm(sortperm(totals, [:agelabel, :total]))

    spectotals = by(samples, :taxon, spectotal = :total =>sum)
    sort!(spectotals, :spectotal)
    topspec = Set(last(spectotals, 5).taxon)

    by(samples, :sample) do df
        df = filter(row-> !in(row.taxon, topspec), df)
        if size(df, 1) == 0
            1
        else
            other = sum(df.total)
            sample = first(df).sample
            push!(samples, (sample=sample, taxon="other", total=other, agelabel=agelabels[sample]))
            1
        end
    end
    push!(topspec, "other")
    filter!(row-> row.taxon in topspec, samples)


    for s in unique(samples.sample)
        missing_tax = setdiff(topspec, samples[samples.sample .== s, :taxon])
        for t in missing_tax
            push!(samples, (sample=s, taxon=t, total=0., agelabel=agelabels[s]))
        end
    end

    xs = dictionary(s=> i for (s,i) in eachrow(totals[!,[:sample,:ord]]))
    samples.x = [xs[s] for s in samples.sample]
    sort!(samples, :x)
    samples
end

gabad = let geneset = Set(nko["GABA degradation"])
    filt = filter(row-> row.ko in geneset, koslong)
    samples = by(filt, [:sample, :taxon]) do sample
        (total=sum(sample.abundance), agelabel=first(sample).agelabel)
    end

    sort!(samples, [:agelabel, :sample])
    totals = filter(row-> isnothing(row.taxon), samples)
    filter!(row-> !isnothing(row.taxon), samples)
    totals.ord = invperm(sortperm(totals, [:agelabel, :total]))

    spectotals = by(samples, :taxon, spectotal = :total =>sum)
    sort!(spectotals, :spectotal)
    topspec = Set(last(spectotals, 5).taxon)

    by(samples, :sample) do df
        df = filter(row-> !in(row.taxon, topspec), df)
        if size(df, 1) == 0
            1
        else
            other = sum(df.total)
            sample = first(df).sample
            push!(samples, (sample=sample, taxon="other", total=other, agelabel=agelabels[sample]))
            1
        end
    end
    push!(topspec, "other")
    filter!(row-> row.taxon in topspec, samples)


    for s in unique(samples.sample)
        missing_tax = setdiff(topspec, samples[samples.sample .== s, :taxon])
        for t in missing_tax
            push!(samples, (sample=s, taxon=t, total=0., agelabel=agelabels[s]))
        end
    end

    xs = dictionary(s=> i for (s,i) in eachrow(totals[!,[:sample,:ord]]))
    samples.x = [xs[s] for s in samples.sample]
    sort!(samples, :x)
    samples
end

gabas = let geneset = Set(nko["GABA synthesis"])
    filt = filter(row-> row.ko in geneset, koslong)
    samples = by(filt, [:sample, :taxon]) do sample
        (total=sum(sample.abundance), agelabel=first(sample).agelabel)
    end

    sort!(samples, [:agelabel, :sample])
    totals = filter(row-> isnothing(row.taxon), samples)
    filter!(row-> !isnothing(row.taxon), samples)
    totals.ord = invperm(sortperm(totals, [:agelabel, :total]))

    spectotals = by(samples, :taxon, spectotal = :total =>sum)
    sort!(spectotals, :spectotal)
    topspec = Set(last(spectotals, 5).taxon)

    by(samples, :sample) do df
        df = filter(row-> !in(row.taxon, topspec), df)
        if size(df, 1) == 0
            1
        else
            other = sum(df.total)
            sample = first(df).sample
            push!(samples, (sample=sample, taxon="other", total=other, agelabel=agelabels[sample]))
            1
        end
    end
    push!(topspec, "other")
    filter!(row-> row.taxon in topspec, samples)


    for s in unique(samples.sample)
        missing_tax = setdiff(topspec, samples[samples.sample .== s, :taxon])
        for t in missing_tax
            push!(samples, (sample=s, taxon=t, total=0., agelabel=agelabels[s]))
        end
    end

    xs = dictionary(s=> i for (s,i) in eachrow(totals[!,[:sample,:ord]]))
    samples.x = [xs[s] for s in samples.sample]
    sort!(samples, :x)
    samples
end
