using Pkg
Pkg.activate("./")
using JLD2
using CSV
using DataFrames
using MultivariateStats
using Statistics
using ECHOAnalysis

@load "analysis/figures/assets/taxa.jld2" speciesmds allmeta speciesmdsaxes
@load "analysis/figures/assets/figure1b.jld2" kidsspeciesmds kidsspeciesmdsaxes allkidsmeta
@load "analysis/figures/assets/figure1c.jld2" unirefaccessorymds ubothmeta unirefaccessorymdsaxes
@load "analysis/figures/assets/figure1d.jld2" kidsunirefaccessorymds kidsunirefaccessorymdsaxes
@load "analysis/figures/assets/figure1e.jld2" r2 r2m qa allpermanovas
@load "analysis/figures/assets/figure1g.jld2" allfsea mdcors
@load "analysis/figures/assets/s1.jld2" speciesdiffs unirefaccessorydiffs kosdiffs pfamsdiffs

# Get samples for kids under ~14 months

babies = filter(row->!ismissing(row.correctedAgeDays) && row.correctedAgeDays / 365 < 1.2, allmeta)
samples = Set(babies.sample)


# Get files for those samples

fastqs = String[]
for (root, dirs, files) in walkdir("/lovelace/echo/analysis/engaging/")
    !occursin(r"batch\d{3}", root) && continue
    !occursin(r"kneaddata", root) && continue
    filter!(f-> occursin(r"kneaddata\.fastq.gz$", f) &&
                sampleid(stoolsample(f)) in samples,
            files)

    append!(fastqs, joinpath.(root,files))
end

open("/home/kevin/Desktop/babies.txt", "w") do io
    for f in fastqs
        write(io, f * '\n')
    end
end


# ## Submitting jobs to cluster
#
# batch_job.sh is:
# ```
# #!/bin/bash
# #SBATCH -n 1                # Number of cores
# #SBATCH -N 1                # Ensure that all cores are on one machine
# #SBATCH -t 0-08:00          # Runtime in D-HH:MM, minimum of 10 minutes
# #SBATCH -p newnodes         # Partition to submit to
# #SBATCH --mem=8000          # Memory pool for all cores (see also --mem-per-cpu)
# #SBATCH -o /nobackup1/vklepacc/echo/metaphlan3/%j.out
# #SBATCH -e /nobackup1/vklepacc/echo/metaphlan3/%j.err
#
# source /home/vklepacc/.bashrc
# source activate metaphlan3
# metaphlan --input_type fastq ${infile} -o ./output/${sid}_profile.tsv --bowtie2out ./output/${sid}_bowtie2.bt
# ```

# Then, on engaging cluster, run
#
# ```julia
# cd("/nobackup1/vklepacc/echo/metaphlan3/")
#
# fastqs = joinpath.("fastqs", readdir("fastqs"))
#
# sids = map(fastqs) do f
#     m = match(r"fastqs\/(C[\w\-]+)_S\d{1,2}_kneaddata\.fastq\.gz", f)
#     isnothing(m) && error("no match for $f")
#     string(m.captures[1])
# end
#
# for (f, i) in zip(fastqs, sids)
#     @async run(`sbatch --export=infile=$f,sid=$i batch_job.sh`)
# end
# ```

# ## Meanwhile, back on `ada`

newprofiles = readdir("/lovelace/echo/analysis/metaphlan3/profiles/", join=true)
df = DataFrame(CSV.File(newprofiles[1], delim='\t', header=["clade", "taxid", "abundance"], select=[1,2,3], skipto=5))

size(df)
df[!,4]
names(df)
dfs = [DataFrame(CSV.File(f, delim='\t')) for f in newprofiles]
