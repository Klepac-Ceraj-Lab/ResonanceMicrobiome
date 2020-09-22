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

function kneaddata(sample, filelist, outpath)
    path, io = mktemp()
    catfile = path * ".fastq.gz"
    run(pipeline(`cat $filelist`, stdout=catfile))
    cmd = `kneaddata --input $catfile --reference-db /babbage/biobakery_databases/kneaddata/hg37 --reference-db /babbage/biobakery_databases/kneaddata/silva --output $outpath --output-prefix $(sample)_kneaddata --trimmomatic /opt/miniconda3/envs/bb3/share/trimmomatic-0.39-1/ --threads 8`
    run(cmd)
    return
end

for row in eachrow(redo[2:end, :])
    s = row.sample
    isfile(joinpath(kneaddir, "$(s)_kneaddata.fastq")) && continue
    
    @info "running kneaddata on $s"
    raw = filter(f-> occursin(replace(s, "_"=> "-"), f), rawfastq_files)
    try 
        kneaddata(s, raw, kneaddir)
    catch e
        @warn "$s threw error" e
    end
end