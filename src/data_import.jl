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

function functional_profiles(kind; stratify=false, commonfilter=true)
    stratify && error("Functions with species stratification doesn't work yet")
    df = Arrow.Table(@datadep_str "$kind/$kind.arrow") |> DataFrame
    !stratify && (df = filter(row-> !occursin("|", row.feature), df))
    # df[!, :gftemp] = [GeneFunction("") for _ in 1:nrow(df)]
    # for row in eachrow(df)
    #     pieces = split(row.feature, "|")
    #     if length(pieces) == 1
    #         gf = GeneFunction(first(pieces))
    #     else
    #         length(pieces) == 2 || throw(ErrorException("Expected gene funcction, got $pieces"))
            
    #         sp = split(pieces[2], ".")
    #         if first(sp) == "unclassified"
    #             gf = GeneFunction(first(pieces))
    #         else            
    #             length(sp) != 2 && error("something went wrong with $sp")
    #             sp = replace(sp[2], "s__" => "")
    #             sp = Taxon(sp, :species)
    #             gf = GeneFunction(pieces[1], sp)
    #         end
    #     end
    #     row.gftemp = gf
    # end
    # select!(df, Not("feature"))
    # rename!(df, :gftemp=>:feature)

    samples = unique(df.sample)
    commonfilter && sample_filter!(samples)
    gfs = GeneFunction.(unique(df.feature))

    gfmap = Dict(t => i for (i, t) in enumerate(gfs))
    sampmap = Dict(s => i for (i, s) in enumerate(samples))
    mat = spzeros(length(gfs), length(samples))
    
    for row in eachrow(df)
        if commonfilter
            sample_filter(row.sample) || continue
        end
        mat[gfmap[GeneFunction(row.feature)], sampmap[row.sample]] = row.abundance ./ 100
    end
    return CommunityProfile(mat, gfs, MicrobiomeSample.(samples))
end
