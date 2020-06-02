include("../scripts/startup_loadpackages.jl")

@load "analysis/figures/assets/metadata.jld2" allmeta ubothmeta ukidsmeta allkidsmeta allmoms allkids umoms ukids oldkids uboth
allkidsmeta.sample = [String(s) for s in allkidsmeta.sample]
@load "analysis/figures/assets/taxa.jld2" species speciesmds speciesmdsaxes ubothspeciesmds ubothspeciesmdsaxes ukidsspeciesmds ukidsspeciesmdsaxes
@load "analysis/figures/assets/unirefs.jld2" unirefaccessorymds unirefaccessorymdsaxes ubothunirefaccessorymds ubothunirefaccessorymdsaxes ukidsunirefaccessorymds ukidsunirefaccessorymdsaxes
@load "analysis/figures/assets/otherfunctions.jld2" kos kosdiffs kosdm ecs ecsdm pfams pfamsdiffs pfamsdm
@load "analysis/figures/assets/permanovas.jld2" r2 r2m qa allpermanovas species_permanovas unirefaccessory_permanovas kos_permanovas pfams_permanovas
@load "analysis/figures/assets/fsea.jld2" allfsea mdcors
@load "analysis/figures/assets/difs.jld2" speciesdiffs unirefaccessorydiffs kosdiffs pfamsdiffs
@load "analysis/figures/assets/stratkos.jld2" stratkos
@load "analysis/figures/assets/cogquartiles.jld2" quartmeta quartspecies quartspeciesdm quartspeciesmds quartspeciesmdsaxes quartiletests

## 

open("analysis/tables/table1.tsv", "w") do io
    println(io, "thing\tvalue")
    println(io, "Subjects (n)\t$(sum(!=("mom"), ubothmeta.ageLabel))")
    println(io, "Under 1yo (n)\t$(sum(==("1 and under"), ubothmeta.ageLabel))")
    println(io, "Over 2yo (n)\t$(sum(==("2 and over"), ubothmeta.ageLabel))")

    println(io, "With highresolution scan (n)\t$(sum(!ismissing, ubothmeta.hires_total))")
    println(io, "With cognitive function score (n)\t$(sum(!ismissing, ubothmeta.cogScore))")
    println(io, "Both scan and cognitive function (n)\t$(sum(row-> all(!ismissing, (row.hires_total, row.cogScore)), eachrow(ubothmeta)))")

    println(io, "Non-white (%)\t$(round(mean(!=("Caucasian / White"), skipmissing(ubothmeta.simple_race)) *100, digits=2))")
    println(io, "Mixed race (%)\t$(round(mean(r-> occursin(r"[Mm]ixed", r), skipmissing(ubothmeta.simple_race)) *100, digits=2))")
    println(io, "Age in years (mean, SD)\t$(round(mean(skipmissing(ukidsmeta.correctedAgeDays ./ 365)), digits=2)), $(round(std(skipmissing(ukidsmeta.correctedAgeDays ./ 365)), digits=2))")
    println(io, "BMI (mean, SD)\t$(round(mean(skipmissing(ukidsmeta.childBMI)), digits=2)), $(round(std(skipmissing(ukidsmeta.childBMI)), digits=2))")
    println(io, "Maternal SES (mean, SD)\t$(round(mean(skipmissing(ukidsmeta.mother_HHS)), digits=2)), $(round(std(skipmissing(ukidsmeta.mother_HHS)), digits=2))")
end


## ## Supplementary Tables

CSV.write("analysis/tables/supptable1_allpermanovas.csv", allpermanovas)
CSV.write("analysis/tables/supptable2_quartiletests.csv", quartiletests)
CSV.write("analysis/tables/supptable3_allfsea.csv", allfsea[!, Not(:cors)])
CSV.write("analysis/tables/supptable4_allmetadata.csv", 
    select(ubothmeta, [
        :sample,
        :subject,
        :timepoint,
        :correctedAgeDays,
        :breastfeeding,
        :breastFedPercent,
        :birthType,
        :childGender,
        :simple_race,
        :childBMI,
        :mother_HHS,
        :cogScore,
        :cogAssessment,
        :hires_total,
        :limbic_normed,
        :subcortical_normed,
        :neocortical_normed,
        :cerebellar_normed,
    ]))

## Other numbers

println("Total stool: ", nrow(ubothmeta))
println("Has breastfeeding: ", count(!ismissing, ukidsmeta.breastfeeding))
println("Has birthType: ", count(!ismissing, ukidsmeta.birthType))
println("Has SES: ", count(!ismissing, ukidsmeta.mother_HHS))
println("Has BMI: ", count(!ismissing, ukidsmeta.childBMI))
println("N top quartile: ", count(==("top 25%"), quartmeta.quartile))
println("N bottom quartile: ", count(!=("top 25%"), quartmeta.quartile))
println("Minimum age (days): ", minimum(skipmissing(ukidsmeta.correctedAgeDays)))
println("Minimum age (days) with MRI: ", minimum(filter(row-> !ismissing(row.correctedAgeDays) && !ismissing(row.hires_total), ukidsmeta).correctedAgeDays))
println("B. longum prevalence ", count(>(0), vec(occurrences(view(species, sites=ukids, species=["Bifidobacterium_longum"])))) / sum(ukids))

println("Has bf < 6 months: ", 
    nrow(filter(ukidsmeta) do row
        !ismissing(row.correctedAgeDays) &&
        !ismissing(row.breastfeeding) &&
        row.correctedAgeDays < 365/2
    end)
    )
println("Has bf < 1 year: ", 
    nrow(filter(ukidsmeta) do row
        !ismissing(row.correctedAgeDays) &&
        !ismissing(row.breastfeeding) &&
        row.correctedAgeDays < 365
    end)
    )