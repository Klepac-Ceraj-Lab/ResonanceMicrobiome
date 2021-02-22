function resonance_metadata(; commonfilter=true)
    samplesdf = CSV.read(datadep"sample_metadata/sample_metadata.csv", DataFrame)
    if commonfilter
        filter!(row-> sample_filter(row.sample))
    end    
    clinicaldf = CSV.read(datadep"clinical_metadata/clinical_metadata.csv", DataFrame)

    return leftjoin(samplesdf, clinicaldf, on=[:subject, :timepoint])
end

function resonance_metadata(samples; commonfilter=true)
    outdf = resonance_metadata(commonfilter=commonfilter)
    samplemap = Dict(s=>i for (i,s) in enumerate(outdf.sample))
    idx = map(samples) do s
        i = get(samplemap, s, nothing)
        isnothing(i) && throw(ArgumentError("$s not found in metadata"))
        i
    end
    return outdf[idx, :]
end


function taxonomic_profiles(level; commonfilter=true)
    df = Arrow.Table(datadep"taxonomic_profiles/taxonomic_profiles.arrow") |> DataFrame
    df.taxon = [last(BiobakeryUtils._split_clades(c)) for c in df.clade]
    keep = level == :all ? Colon() : [ismissing(c) || c == level for c in clade.(df.taxon)]
    df = df[keep, :]
    
    samples = unique(df.sample)
    commonfilter && sample_filter!(samples)
    taxa = unique(df.taxon)

    taxmap = Dict(t => i for (i, t) in enumerate(taxa))
    sampmap = Dict(s => i for (i, s) in enumerate(samples))
    mat = spzeros(length(taxa), length(samples))
    
    for row in eachrow(df)
        if commonfilter
            sample_filter(row.sample) || continue
        end
        mat[taxmap[row.taxon], sampmap[row.sample]] = row.abundance ./ 100
    end
    return CommunityProfile(mat, taxa, MicrobiomeSample.(samples))
end

function functional_profiles(kind)
    df = Arrow.Table(@datadep_str "$kind/$kind.arrow") |> DataFrame
    
end