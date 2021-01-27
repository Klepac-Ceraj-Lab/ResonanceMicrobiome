function resonance_metadata()
    samplesdf = CSV.read(datadep"sample_metadata/sample_metadata.csv", DataFrame)
    clinicaldf = CSV.read(datadep"clinical_metadata/clinical_metadata.csv", DataFrame)
    return leftjoin(samplesdf, clinicaldf, on=[:subject, :timepoint])
end

function resonance_metadata(samples)
    samplesdf = CSV.read(datadep"sample_metadata/sample_metadata.csv", DataFrame)
    clinicaldf = CSV.read(datadep"clinical_metadata/clinical_metadata.csv", DataFrame)
    outdf = leftjoin(samplesdf, clinicaldf, on=[:subject, :timepoint])
    samplemap = Dict(s=>i for (i,s) in enumerate(outdf.sample))
    idx = map(samples) do s
        i = get(samplemap, s, nothing)
        isnothing(i) && throw(ArgumentError("$s not found in metadata"))
        i
    end
    return outdf[idx, :]
end


function taxonomic_profiles(level)
    df = Arrow.Table(datadep"taxonomic_profiles/taxonomic_profiles.arrow") |> DataFrame
    df.taxon = [last(BiobakeryUtils._split_clades(c)) for c in df.clade]
    keep = level == :all ? Colon() : [ismissing(c) || c == level for c in clade.(df.taxon)]
    df = df[keep, :]
    
    samples = unique(df.sample)
    taxa = unique(df.taxon)

    taxmap = Dict(t => i for (i, t) in enumerate(taxa))
    sampmap = Dict(s => i for (i, s) in enumerate(samples))
    mat = spzeros(length(taxa), length(samples))
    
    for row in eachrow(df)
        mat[taxmap[row.taxon], sampmap[row.sample]] = row.abundance ./ 100
    end
    return CommunityProfile(mat, taxa, MicrobiomeSample.(samples))
end