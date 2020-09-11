using ECHOAnalysis
using DataFrames
using CSV
import Pkg.TOML: parsefile 

config = parsefile("Data.toml")
allmeta = CSV.File(config["tables"]["joined_metadata"], pool=false) |> DataFrame
sample2batch = Dict(s => "batch$(lpad(b,3,'0'))" for (s,b) in eachrow(select(allmeta, :sample,:batch)))

taxonomic_profiles = readdir(config["filepaths"]["taxonomic_profiles"])
functional_profiles = filter(p-> occursin("genefamilies_relab",p), readdir(config["filepaths"]["functional_profiles"]))

missingtax = setdiff(allmeta.sample, sampleid.(stoolsample.(taxonomic_profiles)))
missingfunc = setdiff(allmeta.sample, sampleid.(stoolsample.(functional_profiles)))

redo = filter(row-> row.sample in missingfunc || row.sample in missingtax, allmeta)

rawfastq = unique(sampleid.(stoolsample.(filter(r-> !occursin("Zymo", r), readdir(config["filepaths"]["rawfastq"])))))

for row in eachrow(redo)
    s = row.sample
    if !(s in rawfastq)
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
    allraw = readdir(config["filepaths"]["rawfastq"])
    
    raw = filter(f-> occursin(replace(s, "_"=> "-"), f), allraw)
    length(raw) != 8 && @warn "$s has more than 8 raw files"

    
end