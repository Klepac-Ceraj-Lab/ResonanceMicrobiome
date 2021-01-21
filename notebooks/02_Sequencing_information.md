# Sequencing data

Stool samples were collected by parents in OMR-200 tubes,
and DNA extraction was performed at Wellesley College (Wellesley, MA).

Nucleic acids were extracted from stool samples
using the RNeasy PowerMicrobiome kit automated on the QIAcube,
excluding the DNA degradation steps.
Extracted DNA was sequenced at the Integrated Microbiome Resource (IMR, Dalhousie University, NS, Canada).

For shotgun metagenomic sequencin, a pooled library (max 96 samples per run)
was prepared using the Illumina Nextera Flex Kit for MiSeq and NextSeq from 1 ng of each sample.
Samples were then pooled onto a plate and sequenced on the Illumina NextSeq 550 platform
using 150+150 bp paired-end “high output” chemistry,
generating ~400 million raw reads and ~120 Gb of sequence. 

## Sequence Read Archive (SRA) submission

```julia
using ResonanceMicrobiome # rexports CSV and DataFrames

reads = CSV.read(datadep"kneaddata_counts/allcounts.csv", DataFrame)
metadata = resonance_metadata()
```

```julia
samplenames = [replace(s, r"_S\d{1,2}_kneaddata"=>"") for s in reads.Sample]
sra_files = ["$(sample)_pair1.fastq.gz" for sample in samplenames]
```


Traits for SRA submission

```julia
sra = DataFrame(
    "Sample Name" => samplenames,
    "file1"       => sra_files,
    "file2"       => replace.(sra_files, Ref("pair1"=>"pair2")),
)

sra[!, "organism"]            .= "human gut metagenome" # txid408170
sra[!, "host"]                .= "Homo sapiens"
sra[!, "investigation_type"]  .= "metagenome"
sra[!, "project_name"]        .= "Human gut metagenomes of healthy mothers and their children, Jan 19 '21"
sra[!, "lat_lon"]             .= "41°48'37\"N, 71°24'24\"W" # location of RI Women and Infants hospital
sra[!, "geo_loc_name"]        .= "/country=\"USA: Rhode Island\""
sra[!, "env_package"]         .= "human-gut"
sra[!, "specific_host"]       .= "human"
sra[!, "rel_to_oxygen"]       .= "anaerobic"
sra[!, "samp_collect_device"] .= "Omnigene OMR-200 tube"
sra[!, "seq_meth"]            .= "illumina NextSeq 550"
sra[!, "seq_quality_check"]   .= "none"
sra[!, "host_disease_stat"]   .= "none"
sra[!, "host_body_product"]   .= "fma64183"
sra[!, "samp_store_temp"]     .= "-80 degree celcius"
sra[!, "nucl_acid_ext"]       .= "https://www.qiagen.com/us/resources/resourcedetail?id=84c1f2e7-8db6-4957-a504-92bf9f82dd84"
```

And now for sample-specific attributes:

```julia
# make DataFrame with same rows as SRA
samplemeta = DataFrame(sample = replace.(sra."Sample Name", Ref("-"=>"_")))
samplemeta = leftjoin(samplemeta, metadata, on=:sample)

# sra[!, "adapters"]           .=
# sra[!, "mid"]                .= {Multiplex ID}
sra[!, "collection_date"] = samplemeta.DOC
sra[!, "host_age"]        = map(a-> ismissing(a) ? missing : string.(a) .* " days", samplemeta.correctedAgeDays)
sra[!, "host_sex"]        = map(row-> startswith(row.sample, "M") ? "Female" : row.childGender, eachrow(samplemeta))
```

Generate output for BioSample Attributes

```julia
CSV.write("data/biosample_attributes.tsv", sra[!, Not([:file1, :file2])])
```
## Visualizing read data

```julia
using CairoMakie

figure = Figure(resolution = (1200, 800))

f1 = figure[1,1] = Axis(title="Raw vs Final Reads")
```