"""
    airtable_metadata(key=ENV["AIRTABLE_KEY"])
Get fecal sample metadata table from airtable.
The API `key` comes from https://airtable.com/account.

This is unlikely to work if you're not in the VKC lab,
but published sample metadata is available from OSF.io
using `datadep"sample metadata"`.
"""
function airtable_metadata(key=Airtable.Credential())
    records = []
    req = Airtable.get(key, "/v0/appyRaPsZ5RsY4A1h", "Master"; view="ALL_NO_EDIT", filterByFormula="NOT({Mgx_batch}='')")
    append!(records, req.records)
    while haskey(req, :offset) && length(records) < 2200
        @info "Making another request"
        req = Airtable.get(key, "/v0/appyRaPsZ5RsY4A1h/", "Master"; view="ALL_NO_EDIT", filterByFormula="NOT({Mgx_batch}='')", offset=req.offset)
        append!(records, req.records)
        sleep(0.250)
    end

    df = DataFrame()
    for record in records
        append!(df, filter(p -> !(last(p) isa AbstractArray), record.fields), cols=:union)
    end

    rename!(df, "SampleID"=>"sample", "TimePoint"=>"timepoint", "SubjectID"=>"subject")
    
    transform!(df, "subject"   => ByRow(s-> parse(Int, s)) => "subject",
                   "timepoint" => ByRow(tp-> parse(Int, tp)) => "timepoint",
                   "Mgx_batch" => ByRow(b-> !ismissing(b) ? parse(Int, match(r"Batch (\d+)", b).captures[1]) : missing) => "Mgx_batch",
                   "16S_batch" => ByRow(b-> !ismissing(b) ? parse(Int, match(r"Batch (\d+)", b).captures[1]) : missing) => "16S_batch")
    return select(df, Cols(:sample, :subject, :timepoint, :))
end

"""
    sample_filter(sample(s))

Set of common filters to remove samples from analysis.
Normalizes sample names by replacing "-" with "_".

- For individual items, returns true if it should be kept.
- For arrays, returns a filtered array (see also [`sample_filter!`](@ref))
- For community progfiles, returns a filtered community profile.

**Currently removes:**

- anything that doesn't match `r"[CM]\\d{4}_\\d+F_1A`.
  This incidentally eleminates:
  - ethanol samples (which have `r"_\\d+E_"' in their name)
  - anything that's not the first replicate / aliquot (should have `_1A` at the end)
- the samples `("M1295_3F_1A", "M1322_2F_1A", "C1155_4F_1A", "M1367_2F_1A")`,
  which failed QC on the sequencing side
"""
function sample_filter(sample::AbstractString)
    sample = replace(sample, "-"=>"_")
    !occursin(r"[CM]\d{4}_\d+F_1A", sample) && return false
    in(sample, ("M1295_3F_1A", "M1322_2F_1A", "C1155_4F_1A", "M1367_2F_1A")) && return false
    return true
end

sample_filter(sample::Symbol) = sample_filter(String(sample))
sample_filter(sample::AbstractSample) = sample_filter(name(sample))

sample_filter(samples::AbstractArray) = filter(sample_filter, samples)
sample_filter(comm::CommunityProfile) = comm[:, sample_filter(samplenames(comm))]

"""
    sample_filter!(samples::AbstractArray)

Mutating version of [`sample_filter`](@ref)
"""
sample_filter!(samples::AbstractArray) = filter!(sample_filter, samples)

function post_fetch_clinical(csv)
    clinical = CSV.read(csv, DataFrame)
    rename!(clinical, "subjectID"=> "subject")
    filter!(row-> all(!ismissing, [row.subject, row.timepoint]), clinical)
    disallowmissing!(clinical, [:subject, :timepoint])
    CSV.write("clinical_metadata.csv", clinical)
end


function post_fetch_knead(tarball)
    unpack(tarball)
    allcounts = DataFrame()
    for (root, dir, files) in walkdir(pwd())
        length(files) == 1 || continue
        counts = first(files)
        counts == "kneaddata_read_counts.txt" || continue
        batch = match(r"biobakery3\/batch(\d{3})\/output", root).captures[1]
        batch = parse(Int, batch)
        df = CSV.read(joinpath(root, counts), DataFrame)
        df[!, :batch] .= batch
        filter!(row-> sample_filter(row.Sample), df)
        append!(allcounts, df)
    end
    CSV.write("allcounts.csv", allcounts)
    rm("biobakery3", recursive=true, force=true)
    rm("download", force=true)  # we don't know where this file comes from but we want it gone
end

function post_fetch_taxa(tarball)
    unpack(tarball)
    profile_paths = String[]
    for (root, dir, files) in walkdir(pwd())
        filter!(f-> occursin("profile.tsv", f), files)
        append!(profile_paths, joinpath.(root, files))
    end
    
    df = DataFrame()
    for profile in profile_paths
        pdf = CSV.read(profile, DataFrame, datarow=5, header=["clade", "NCBI_taxid", "abundance", "additional_species"])
        pdf[!,:sample] .= replace(replace(first(splitext(basename(profile))), r"_S\d{1,2}_profile"=> ""), "-"=> "_")
        append!(df, pdf)
    end
    Arrow.write("taxonomic_profiles.arrow", df)

    # profiles = Tables.partitioner(profile_paths) do profile
    #     df = CSV.File(profile, datarow=5, header=["clade", "NCBI_taxid", "abundance", "additional_species"])
    # end
    # Arrow.write("taxonomic_profiles.arrow", profiles)
    rm("biobakery3", recursive=true, force=true)
    rm("download", force=true)  # we don't know where this file comes from but we want it gone
end

function _process_osf_functional_profiles()
    profile_paths = String[]
    for (root, dir, files) in walkdir(pwd())
        filter!(f-> occursin(".tsv", f), files)
        append!(profile_paths, joinpath.(root, files))
    end
    
    df = DataFrame()
    for profile in profile_paths
        pdf = CSV.read(profile, DataFrame, header=["feature", "abundance"], datarow=2)
        pdf[!,:sample] .= replace(replace(first(splitext(basename(profile))), r"_S\d{1,2}_.+"=> ""), "-"=> "_")
        transform(pdf, :feature => ByRow(s-> begin
            things = split(s, "|")
            (p1, p2) = length(things) == 1 ? (first(things), nothing) : things
            tax = isnothing(p2) ? missing : replace(p2, r"g__\w+\.s__"=>"")
            (p1, tax)
        end) => [:feature, :species])
        
        append!(df, pdf)
    end
    return df
end


function post_fetch_kos(tarball)
    unpack(tarball)
    
    df = _process_osf_functional_profiles()
    Arrow.write("kos.arrow", df)

    # profiles = Tables.partitioner(profile_paths) do profile
    #     df = CSV.File(profile, datarow=5, header=["clade", "NCBI_taxid", "abundance", "additional_species"])
    # end
    # Arrow.write("taxonomic_profiles.arrow", profiles)
    rm("biobakery3", recursive=true, force=true)
    rm("download", force=true)  # we don't know where this file comes from but we want it gone
end


function post_fetch_unirefs(tarball)
    unpack(tarball)
    
    df = _process_osf_functional_profiles()
    Arrow.write("unirefs.arrow", df)

    # profiles = Tables.partitioner(profile_paths) do profile
    #     df = CSV.File(profile, datarow=5, header=["clade", "NCBI_taxid", "abundance", "additional_species"])
    # end
    # Arrow.write("taxonomic_profiles.arrow", profiles)
    rm("biobakery3", recursive=true, force=true)
    rm("download", force=true)  # we don't know where this file comes from but we want it gone
end

function post_fetch_pfams(tarball)
    unpack(tarball)
    
    df = _process_osf_functional_profiles()
    Arrow.write("pfams.arrow", df)

    # profiles = Tables.partitioner(profile_paths) do profile
    #     df = CSV.File(profile, datarow=5, header=["clade", "NCBI_taxid", "abundance", "additional_species"])
    # end
    # Arrow.write("taxonomic_profiles.arrow", profiles)
    rm("biobakery3", recursive=true, force=true)
    rm("download", force=true)  # we don't know where this file comes from but we want it gone
end

function post_fetch_ecs(tarball)
    unpack(tarball)
    
    df = _process_osf_functional_profiles()
    Arrow.write("ecs.arrow", df)

    # profiles = Tables.partitioner(profile_paths) do profile
    #     df = CSV.File(profile, datarow=5, header=["clade", "NCBI_taxid", "abundance", "additional_species"])
    # end
    # Arrow.write("taxonomic_profiles.arrow", profiles)
    rm("biobakery3", recursive=true, force=true)
    rm("download", force=true)  # we don't know where this file comes from but we want it gone
end


function post_fetch_uniprot(tarball)
    unpack(tarball)
    rm("download", force=true)
end

post_fetch_osf(::Any) = rm("download", force=true)
post_fetch_osf() = rm("download", force=true)

function __init__()
    register(DataDep(
        "taxonomic_profiles",
        """
        Taxonomic profiles generated by MetaPhlAn.
        """,
        "https://osf.io/kyh48/download",
        "35a172c44e41ea4a7aaa51ce3e94526accc5a7e7f03e44e6c553bd8f0915ee87";
        post_fetch_method = ResonanceMicrobiome.post_fetch_taxa
        )
    )

    register(DataDep(
        "kos",
        """
        Functional profiles generated by humann.
        """,
        "https://osf.io/736sv/download",
        "2d01ecb311b72381e50d32d8082f4cb5b6aa343688cef2793af6ef4f770571f5";
        post_fetch_method = ResonanceMicrobiome.post_fetch_kos
        )
    )
    
    register(DataDep(
        "unirefs",
        """
        Functional profiles generated by humann.
        """,
        "https://osf.io/wfbxr/download",
        "53b72556a15b21e7982beaeb3e262cd43b18bfad013027332181330c7c169ce1";
        post_fetch_method = ResonanceMicrobiome.post_fetch_unirefs
        )
    )

    register(DataDep(
        "pfams",
        """
        Functional profiles generated by humann.
        """,
        "https://osf.io/eu48m/download",
        "2339dc1961adb65af7c8219656e7d81d36419ea0e097e077c1a728fcfc7dda54";
        post_fetch_method = ResonanceMicrobiome.post_fetch_pfams
        )
    )

    register(DataDep(
        "ecs",
        """
        Functional profiles generated by humann.
        """,
        "https://osf.io/jncwk/download",
        "9bec037b9f60a687bef3e497b6ce1a1942c3500153f5174adf25818be264f0dd";
        post_fetch_method = ResonanceMicrobiome.post_fetch_ecs
        )
    )

    register(DataDep(
        "kneaddata_counts",
        """
        Diagnostic information about sequence reads after cleaning with KneadData.
        """,
        "https://osf.io/ugxjs/download",
        "ad75c18b4987147c1a282a662195e0b717140ee3aacbab7f412c4c8e15c0a1ad";
        post_fetch_method = ResonanceMicrobiome.post_fetch_knead,
        )
    )


    register(DataDep(
        "uniprot",
        """
        Downloads from uniprot
        """,
        "https://osf.io/z6w95/download",
        "c14420f1dd5a0236f78fcfefc4f21bbb38bdae824d7688b8223d66518427577f";
        post_fetch_method = ResonanceMicrobiome.post_fetch_uniprot,
        )
    )


    register(DataDep(
        "clinical_metadata",
        """
        Clinical and subject-specific metadata.
        """,
        "https://osf.io/53b6p/download",
        "08bb0c77c17e54e08655a3c2c831cd0f1b85a5401624c77cdc75c3a0f64cb930";
        post_fetch_method = ResonanceMicrobiome.post_fetch_clinical
        )
    )
    
    register(DataDep(
        "sample_metadata",
        """
        Sample specific metadata.
        """,
        "https://osf.io/5z8tw/download",
        "a070fb88f7455e63e5cc81f0e8c065dafa2935d7ca6a3cc15e099cd7b923dead";
        post_fetch_method = ResonanceMicrobiome.post_fetch_osf
        )
    )

    register(DataDep(
        "amplicon_list",
        """
        List of file names from V4V5 amplicon sequencing
        """,
        "https://osf.io/zca69/download",
        "27844e4e96325f98ffa583768b256b7baed3d5530d7a8e50e204fdb80e61d50b";
        post_fetch_method = ResonanceMicrobiome.post_fetch_osf

    ))
end