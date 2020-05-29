include("../scripts/startup_loadpackages.jl")

@load "analysis/figures/assets/metadata.jld2" allmeta ubothmeta ukidsmeta allkidsmeta allmoms allkids umoms ukids oldkids uboth
allkidsmeta.sample = [String(s) for s in allkidsmeta.sample]
@load "analysis/figures/assets/taxa.jld2" species speciesmds kidsspeciesmds kidsspeciesmdsaxes
@load "analysis/figures/assets/unirefs.jld2" unirefaccessorymds unirefaccessorymdsaxes kidsunirefaccessorymds kidsunirefaccessorymdsaxes
@load "analysis/figures/assets/otherfunctions.jld2" kos kosdiffs kosdm ecs ecsdm pfams pfamsdiffs pfamsdm
@load "analysis/figures/assets/permanovas.jld2" r2 r2m qa allpermanovas species_permanovas unirefaccessory_permanovas kos_permanovas pfams_permanovas
@load "analysis/figures/assets/fsea.jld2" allfsea mdcors
@load "analysis/figures/assets/difs.jld2" speciesdiffs unirefaccessorydiffs kosdiffs pfamsdiffs
@load "analysis/figures/assets/stratkos.jld2" stratkos
@load "analysis/figures/assets/cogquartiles.jld2" quartmeta quartspecies quartspeciesdm quartspeciesmds quartspeciesmdsaxes #quartiletests
allfsea.median = map(median, allfsea.cors)

## 

kids = view(allmeta, allmeta.ageLabel .!= "mom", :)
ukids = let samples = Set(sampleid.(uniquetimepoints(allmeta.sample, takefirst=true, samplefilter=iskid)))
    map(row-> in(row.sample, samples), eachrow(allmeta))
end
ukidsmeta = view(allmeta, ukids, :)

open("analysis/tables/table1.tsv", "w") do io
    println(io, "thing\tvalue")
    println(io, "All samples (n)\t$(size(allmeta, 1))")
    println(io, "Total subjects (n)\t$(size(ubothmeta, 1))")
    println(io, "Moms (n)\t$(sum(==("mom"), ubothmeta.ageLabel))")
    println(io, "Kids (n)\t$(sum(!=("mom"), ubothmeta.ageLabel))")
    println(io, "Kids under 1yo (n)\t$(sum(==("1 and under"), ubothmeta.ageLabel))")
    println(io, "Kids over 2yo (n)\t$(sum(==("2 and over"), ubothmeta.ageLabel))")

    println(io, "Kids with highres (n)\t$(sum(!ismissing, ubothmeta.hires_total))")
    println(io, "Kids with cogscore (n)\t$(sum(!ismissing, ubothmeta.cogScore))")
    println(io, "Kids with both (n)\t$(sum(row-> all(!ismissing, (row.hires_total, row.cogScore)), eachrow(ubothmeta)))")

    println(io, "Non-white kids (%)\t$(round(mean(!=("Caucasian / White"), skipmissing(ubothmeta.simple_race)) *100, digits=2))")
    println(io, "Mixed race kids (%)\t$(round(mean(r-> occursin(r"[Mm]ixed", r), skipmissing(ubothmeta.simple_race)) *100, digits=2))")
    println(io, "Age in years (mean, SD)\t$(round(mean(skipmissing(ukidsmeta.correctedAgeDays ./ 365)), digits=2)), $(round(std(skipmissing(ukidsmeta.correctedAgeDays ./ 365)), digits=2))")
    println(io, "BMI (mean, SD)\t$(round(mean(skipmissing(ukidsmeta.childBMI)), digits=2)), $(round(std(skipmissing(ukidsmeta.childBMI)), digits=2))")
    println(io, "Maternal SES (mean, SD)\t$(round(mean(skipmissing(ukidsmeta.mother_HHS)), digits=2)), $(round(std(skipmissing(ukidsmeta.mother_HHS)), digits=2))")
end


# ## Supplementary Tables

CSV.write("analysis/tables/supptable1_allpermanovas.csv", allpermanovas)
CSV.write("analysis/tables/supptable2_allfsea.csv", allfsea[!, Not(:cors)])
