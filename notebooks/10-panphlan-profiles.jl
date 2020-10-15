using ECHOAnalysis
using DataFrames
using CSV
using TOML: parsefile
using Conda

ENV["PATH"] = "/home/kevin/miniconda3/bin:" * ENV["PATH"]
config = parsefile("Data.toml")

allmeta = CSV.File(config["tables"]["joined_metadata"], pool=false) |> DataFrame
filter!(:correctedAgeDays=> !ismissing, allmeta)
let samples = Set(sampleid.(uniquetimepoints(stoolsample.(allmeta.sample), takefirst=false, samplefilter=iskid)))
    filter!(:sample => s-> s ∈ samples, allmeta)  
end


sigbugs = CSV.File(joinpath(config["output"]["tables"], "oldkidsquartiletests.csv")) |> DataFrame
sigbugs = sigbugs[sigbugs.qvalue .< 0.2, :species]

s = first(allmeta.sample)
b = first(allmeta.batch)

for bug in sigbugs
    isdir("/babbage/panphlan_pangenomes/$bug") || continue
    for row in eachrow(allmeta)
        s = row.sample
        b = row.batch
        kp = joinpath("/lovelace/echo/analysis/biobakery3/", "batch$(lpad(b, 3, "0"))", "output", "kneaddata")
        kneaded = joinpath.(Ref(kp), filter(f-> occursin(replace(s, "_"=>"-"), f) && occursin("kneaddata_paired", f), readdir(kp)))

        run(`panphlan_map.py -i $(kneaded[1]) -i $(kneaded[2]) -p /babbage/panphlan_pangenomes/$(bug)/$(bug)_pangenome.tsv --indexes /babbage/panphlan_pangenomes/$bug/$bug -o $(joinpath(panphlan_out, bug, "$s.tsv"))`)
    end
end