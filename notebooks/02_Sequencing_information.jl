# # Sequencing data

# Stool samples were collected by parents in OMR-200 tubes,
# and DNA extraction was performed at Wellesley College (Wellesley, MA).

# Nucleic acids were extracted from stool samples
# using the RNeasy PowerMicrobiome kit automated on the QIAcube,
# excluding the DNA degradation steps.
# Extracted DNA was sequenced at the Integrated Microbiome Resource (IMR, Dalhousie University, NS, Canada).

# For shotgun metagenomic sequencin, a pooled library (max 96 samples per run)
# was prepared using the Illumina Nextera Flex Kit for MiSeq and NextSeq from 1 ng of each sample.
# Samples were then pooled onto a plate and sequenced on the Illumina NextSeq 550 platform
# using 150+150 bp paired-end “high output” chemistry,
# generating ~400 million raw reads and ~120 Gb of sequence. 

# ## Sequence Read Archive (SRA) submission

using ResonanceMicrobiome # rexports CSV and DataFrames

metadata = resonance_metadata() # depends on datadep"clinical_metadata" and datadep"sample_metadata"

samples = [replace(s, r"_S\d{1,2}_kneaddata"=>"") for s in metadata.sample]
## remove ethanol samples
filter!(s-> !occursin(r"_\d+E_", s), samples)
## remove replicates
filter!(s-> !occursin(r"_\d+F_2", s), samples)
filter!(s-> !occursin(r"_\d+F_\d[^A]", s), samples)

# remove samples that failed QC
filter!(s-> !in(s, ("M1295_3F_1A", "M1322_2F_1A", "C1155_4F_1A", "M1367_2F_1A")), samples)

sra_files = ["$(replace(sample, "_"=>"-"))_pair1.fastq.gz" for sample in samples]
findall(==("C0016-3F-1A_pair1.fastq.gz"), sra_files)
biosample[findall(==("C0066_6F_1A"), biosample."Sample Name"), :file1]


# Global traits for SRA submission

biosample = DataFrame(
    "Sample Name" => samples,
    "file1"       => sra_files,
    "file2"       => replace.(sra_files, Ref("pair1"=>"pair2")),
)

biosample[!, "organism"]            .= "human gut metagenome" # txid408170
biosample[!, "host"]                .= "Homo sapiens"
biosample[!, "investigation_type"]  .= "metagenome"
biosample[!, "project_name"]        .= "Human gut metagenomes of healthy mothers and their children, Jan 19 '21"
biosample[!, "lat_lon"]             .= "41.8786 N 71.3831 W" # location of Pawtucket, RI
biosample[!, "geo_loc_name"]        .= "USA: Rhode Island"
biosample[!, "env_biome"]           .= "not aplicable"
biosample[!, "env_feature"]         .= "not aplicable"
biosample[!, "env_material"]        .= "not aplicable"
biosample[!, "env_package"]         .= "human-gut"
biosample[!, "rel_to_oxygen"]       .= "anaerobe"
biosample[!, "samp_collect_device"] .= "Omnigene OMR-200 tube"
biosample[!, "seq_meth"]            .= "illumina NextSeq 550"
biosample[!, "seq_quality_check"]   .= "none"
biosample[!, "host_disease_stat"]   .= "none"
biosample[!, "host_body_product"]   .= "fma64183"
biosample[!, "samp_store_temp"]     .= "-80 degree celcius"
biosample[!, "nucl_acid_ext"]       .= "https://www.qiagen.com/us/resources/resourcedetail?id=84c1f2e7-8db6-4957-a504-92bf9f82dd84"

# And now for sample-specific attributes:

## make DataFrame with same rows as SRA
samplemeta = DataFrame(sample = replace.(biosample."Sample Name", Ref("-"=>"_")))
samplemeta = leftjoin(samplemeta, metadata, on=:sample)

## biosample[!, "adapters"]           .=
## biosample[!, "mid"]                .= {Multiplex ID}
biosample[!, "subject_id"] = samplemeta.subject
biosample[!, "timepoint"] = samplemeta.timepoint
biosample[!, "collection_date"] = samplemeta.DOC
biosample[!, "host_age"]        = map(a-> ismissing(a) ? missing : string.(a) .* " days", samplemeta.correctedAgeDays)
biosample[!, "host_sex"]        = map(row-> startswith(row.sample, "M") ? "Female" : row.childGender, eachrow(samplemeta))

# Generate output for BioSample Attributes.
CSV.write("output/02_biosample_attributes.tsv", biosample[!, Not([:file1, :file2])], delim='\t')

sra = select(biosample, ["Sample Name", "file1", "file2"])
findall(==("C0016-3F-2A_pair1.fastq.gz"), sra."file1")
findall(==("C0016-3F-2A_pair1.fastq.gz"), biosample."file1")


sra.title = map(row-> "Shotgun metagenomic sequence of stool sample: subject $(row.subject_id), timepoint $(row.timepoint)", eachrow(biosample))
sra.library_ID = map(s-> "$s-mgx", sra."Sample Name")
sra[!, "library_strategy"] .= "WGS"
sra[!, "library_source"] .= "METAGENOMIC"
sra[!, "library_selection"] .= "RANDOM"
sra[!, "library_layout"] .= "paired"
sra[!, "platform"] .= "ILLUMINA"
sra[!, "instrument_model"] .= "NextSeq 550"
sra[!, "design_description"] .= replace("""
    Libraries were prepared with Illumina Nextera Flex Kit
    for MiSeq and NextSeq from 1 ng of each sample.
    Samples were then pooled onto a plate and sequenced
    on the Illumina NextSeq 550 platform
    using 150+150 bp paired-end “high output” chemistry""",
    '\n'=> " ")
sra[!, "filetype"] .= "fastq"
rename!(sra, "file1"=> "filename", "file2"=>"filename2", "Sample Name"=>"sample_name")

# Generate output for SRA Attributes.

CSV.write("output/02_sra_attributes.tsv", sra[!, :], delim='\t')


# ## Visualizing read data

using CairoMakie
using AbstractPlotting.ColorSchemes
using Chain

colormap = ColorSchemes.seaborn_bright6.colors
reads = CSV.read(datadep"kneaddata_counts/allcounts.csv", DataFrame)
sort!(reads, ["batch", "raw pair1"])

select!(reads, "Sample"=>"sample", "batch"=>"batch", 
                names(reads, r"raw") => (+) => "total",
                names(reads, r"trimmed") => (+) => "trimmed",
                names(reads, r"decon") => (+) => "human",
                names(reads, r"final") => (+) => "final",
                )

## adjust to differences
reads.human = reads.trimmed .- reads.human
reads.trimmed = reads.total .- reads.trimmed

## convert to fractions
reads.human_frac = reads.human ./ reads.total
reads.trimmed_frac = reads.trimmed ./ reads.total
reads.final_frac = reads.final ./ reads.total

figure = Figure(resolution = (1200, 800))

fig_a = figure[1:2,1:2] = Axis(figure, title="Raw vs Final Reads", xlabel="Batch #", ylabel="Reads (n)")
batch_edges = map(unique(reads.batch)) do b
    inds = findall(==(b), reads.batch)
    center = extrema(inds)
end

#-

fig_a.xticks = (batch_ticks, string.(1:14))
tightlimits!(fig_a, Bottom())

bar1 = barplot!(fig_a, 1:nrow(reads), reads."total", color=colormap[1], label="total")
bar2 = barplot!(fig_a, 1:nrow(reads), reads."final", color=colormap[5], label="final")
fig_a_leg = figure[1,1] = Legend(figure, fig_a, halign=:right, valign=:top, margin = (10,10,10,10))

fig_b = figure[1,3] = Axis(figure, title="Read Stats", ylabel="Reads (n)", xticks=(1:3, ["Final", "Trimmed", "human"]))
violin!(fig_b, repeat([1,2,3], inner=nrow(reads)), [reads.final; reads.trimmed; reads.human], color=:gray50)
fig_c = figure[2,3] = Axis(figure, title="Read Stats", ylabel="Reads (n / total)", xticks=(1:3, ["Final", "Trimmed", "human"]))
violin!(fig_c, repeat([1,2,3], inner=nrow(reads)), [reads.final_frac; reads.trimmed_frac; reads.human_frac], color=:gray50)

#-

isdir("figures") || mkdir("figures") # make directory if it doesn't exist
CairoMakie.save("figures/02_read_qc.svg", figure)
