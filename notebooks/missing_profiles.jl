using ECHOAnalysis
using DataFrames
using CSV
import TOML: parsefile 

ENV["PATH"] = "/opt/miniconda3/envs/bb3/bin:" * ENV["PATH"]
outpath = "/babbage/echo/bb3_missing"
outfastq = joinpath(outpath, "rawfastq")
isdir(outpath) || mkdir(outpath)
isdir(outfastq) || mkdir(outfastq)

config = parsefile("Data.toml")
allmeta = CSV.File(config["tables"]["joined_metadata"], pool=false) |> DataFrame
sample2batch = Dict(s => "batch$(lpad(b,3,'0'))" for (s,b) in eachrow(select(allmeta, :sample,:batch)))

taxonomic_profiles = readdir(config["filepaths"]["taxonomic_profiles"])
functional_profiles = filter(p-> occursin("genefamilies_relab",p), readdir(config["filepaths"]["functional_profiles"]))

missingtax = setdiff(allmeta.sample, sampleid.(stoolsample.(taxonomic_profiles)))
missingfunc = setdiff(allmeta.sample, sampleid.(stoolsample.(functional_profiles)))

redo = filter(row-> row.sample in missingfunc || row.sample in missingtax, allmeta)

rawfastq_samples = unique(sampleid.(stoolsample.(filter(r-> !occursin("Zymo", r), readdir(config["filepaths"]["rawfastq"])))))
rawfastq_files = readdir(config["filepaths"]["rawfastq"], join=true)

for row in eachrow(redo)
    s = row.sample
    if !(s in rawfastq_samples)
        @warn "Sample hasn't been sequenced" row.sample
        continue
    end

    if s in missingtax
        b = sample2batch[s]
        b == "batch014" && continue 
        if any(p-> occursin(s, p), readdir(joinpath(config["filepaths"]["biobakery"], b, "output/metaphlan")))
            @info "$s in $b has a taxonomic profile, but is not in $(config["filepaths"]["taxonomic_profiles"])"
        end
    end

    if s in missingfunc
        b = sample2batch[s]
        if any(p-> occursin(s, p), readdir(joinpath(config["filepaths"]["biobakery"], b, "output", "humann", "main")))
            @info "$s in $b has a functional profile, but is not in $(config["filepaths"]["functional_profiles"])"
        end
    end

    
    raw = filter(f-> occursin(replace(s, "_"=> "-"), f), rawfastq_files)
    length(raw) != 8 && @warn "$s has more than 8 raw files"
    @show 
    for f in raw
        o = joinpath(outfastq,basename(f))
        isfile(o) || cp(f, joinpath(outfastq,basename(f)))
    end
    
end

kneaddir = joinpath(outpath, "kneaddata"); isdir(kneaddir) || mkdir(kneaddir)
metaphlandir = joinpath(outpath, "metaphlan"); isdir(metaphlandir) || mkdir(metaphlandir)
humanndir = joinpath(outpath, "humann"); isdir(humanndir) || mkdir(humanndir)kneaddir

rawfastq_files = readdir(outfastq, join=true)
redo_sample = unique(sampleid.(stoolsample.(basename.(rawfastq_files))))
filter!(row-> row.sample in redo_sample, redo)

function kneaddata(sample, filelist, outpath)
    path, io = mktemp()
    @info "concatenating"
    catfile = path * ".fastq.gz"
    run(pipeline(`cat $filelist`, stdout=catfile))
    @info "kneading"
    cmd = `kneaddata --input $catfile --reference-db /babbage/biobakery_databases/kneaddata/hg37 --reference-db /babbage/biobakery_databases/kneaddata/silva --output $outpath --output-prefix $(sample)_kneaddata --trimmomatic /opt/miniconda3/envs/bb3/share/trimmomatic-0.39-1/ --threads 8`
    run(cmd)
    return
end

function metaphlan(sample, outpath)
    kneadin = joinpath(kneaddir, "$(sample)_kneaddata.fastq")
    profile_out = joinpath(outpath, sample*"_profile.tsv")
    bowtie_out = joinpath(outpath, sample*"_bowtie2.tsv")
    sam_out = joinpath(outpath, sample*".sam")
    db = "/babbage/biobakery_databases/metaphlan"

    cmd = `metaphlan $kneadin $profile_out --bowtie2out $bowtie_out --samout $sam_out --input_type fastq --nproc 8 --bowtie2db $db --index mpa_v30_CHOCOPhlAn_201901`
    run(cmd)
    return
end

function humann(sample, outpath)
    kneadin = joinpath(kneaddir, "$(sample)_kneaddata.fastq")
    taxin = joinpath(metaphlandir, "$(sample)_profile.tsv")
    profile_out = joinpath(humanndir, "main")
    db = "/babbage/biobakery_databases/metaphlan"
    cmd = `humann --input $kneadin --taxonomic-profile $taxin --output $profile_out --threads 8 --remove-temp-output --search-mode uniref90 --output-basename $sample --metaphlan-options '--bowtie2db $db --index mpa_v30_CHOCOPhlAn_201901'`
    run(cmd)
end

for row in eachrow(redo)
    s = row.sample
    raw = filter(f-> occursin(replace(s, "_"=> "-"), f), rawfastq_files)
    s == "M0753_1F_1A" && continue
    try
        if !isfile(joinpath(kneaddir, "$(s)_kneaddata.fastq"))
            @info "running kneaddata on $s"
            kneaddata(s, raw, kneaddir)
        end
    catch e
        @warn "$s threw error" e
        continue
    end

    try
        if !isfile(joinpath(metaphlandir, "$(s)_profile.tsv"))
            @info "running metaphlan on $s"
            metaphlan(s, metaphlandir)
        end
    catch e
        @warn "$s threw error" e
        continue
    end

    try
        if !isfile(joinpath(humanndir, "$(s)_genefamilies.tsv"))
            @info "running humann on $s"
            humann(s, humanndir)
        end
    catch e
        @warn "$s threw error" e
        break
    end
end
