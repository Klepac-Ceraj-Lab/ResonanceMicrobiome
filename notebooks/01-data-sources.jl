# ---
# title: "Notebook 1: Data Sources"
# author: "Kevin Bonham, PhD"
# options:
#     line_width : 120
#     wrap : false
# ---

# ## Description
#
# Raw data for this paper can be [found here](https://zenodo.org/record/3633793).
#
# In order to perform the analyses outlined in the notebooks in this project,
# you'll need to download and unpack these files
# and put it into directories that are accessible.
#
# By default, data generated from the `bioBakery` toolset
# (those found in the `batchXXX_analysis_noknead.tar.gz` archives)
# are expected to be placed in `data/engaging/` folders.
# Paths for other data may be changed in the `data.toml` file
# found in the root directory of this repository.
#
# Some of these paths also set where files will be created.

# ## Metadata
#
# Patient metadata is kept in a FileMaker Pro database,
# and export by John Rogers.
#
# If you have downloaded data from Zenodo,
# you may skip ahead to the next notebook.
# These next sections document how subject metadata
# is extracted from the internal FilemakerPro database.

# ## Metagenome data
#
# Compressed quality-scored sequencing files (`.fastq.gz`)
# from the sequencing facility were concatenated
# and run through the bioBakery metagenome workflow.
#
# This repository contains samples from sequencing batches 001-013.
# Metadata about fecal samples collected can be exported from airtable
# if you have an API key, though the version for this paper
# Can be found in the zenodo archive above.

# ### Running the snakemake workflow
#
# [Github repository link.](https://github.com/Klepac-Ceraj-Lab/snakemake_workflows)
#
# A custom snakemake workflow was used on all samples
# with the following software versions:
#
# - kneaddata v0.7.1
# - metaphlan v3.0-beta
# - humann v3.0.0-beta
#
# The following command was run on the `engaging` compute cluster (`eofe7.mit.edu`)
# from MIT (a Centos7 environment).
#
# ```
# $ snakemake -s /home/vklepacc/software/repos/snakemake_workflows/biobakery_all.snakefile \
#     --configfile config.yaml --cluster-config cluster.yaml \
#     --cluster "sbatch -n {cluster.processors} -N {cluster.nodes} -t {cluster.time} --mem {cluster.memory} -o output/logs/{rule}-%j.out -e output/logs/{rule}-%j.err -p newnodes --exclude=node119" \
#     --latency-wait 15
# ```

# ## Analyzed profiles
#
# All taxonomic and functional profiles are expected to be found
# in a local directory, and should be specified in the Data.toml file.

# ## Conversion to Arrow format

using DrWatson
@quickactivate "ResonancePaper"

using CSV
using DataFrames
using Arrow
using TOML
using ProgressMeter
using ECHOAnalysis
using BiobakeryUtils

config = TOML.parsefile("Data.toml")
taxpath = config["filepaths"]["taxonomic_profiles"]
funcpath = config["filepaths"]["functional_profiles"]

filepaths = filter(isstoolsample∘basename, readdir(taxpath, join = true))

Arrow.write("data/test.arrow", Tables.partitioner(filepaths) do file
    @info file
    sample = stoolsample(basename(file))
    tax = CSV.File(file, header=[:taxon, :taxid, :abundance, :additional_species],
                    skipto=5) |> DataFrame
    tax.taxon = map(last ∘ parsetaxa, tax.taxon)
    tax[!, :sample] .= sampleid(sample)
    return select(tax, [:sample, :taxon, :taxid, :abundance])
end)

let file = first(filepaths)
    sample = stoolsample(basename(file))
    tax = CSV.File(file, header=[:taxon, :taxid, :abundance, :additional_species],
                    skipto=5) |> DataFrame
    @info tax.taxon
    transform!(tax, :taxon => ByRow(first ∘ parsetaxa) => :taxonlevel)

    tax[!, :sample] .= sampleid(sample)
    Arrow.write("data/test.arrow", tax)

end

df = DataFrame(a=[("thing", :thing), ("otherthing", :other)])

Arrow.write("data/test.arrow",df)