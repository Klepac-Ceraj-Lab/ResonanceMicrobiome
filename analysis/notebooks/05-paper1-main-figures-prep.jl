# # Main Figures
#
# Start after notebook 2.
#
# ## Setup
#
# This step takes about 20 min

include("../scripts/startup_loadall.jl")
## Figure 1B
speciesdm = pairwise(BrayCurtis(), species)
speciesmds = fit(MDS, speciesdm, distances=true)
speciesmdsaxes = [v / sum(eigvals(speciesmds)) for v in eigvals(speciesmds)]

ukidsspeciesmds = fit(MDS, speciesdm[ukids,ukids], distances=true)
ukidsspeciesmdsaxes = [v / sum(eigvals(ukidsspeciesmds)) for v in eigvals(ukidsspeciesmds)]

#### Figure 1C

unirefaccessorydm = pairwise(BrayCurtis(), unirefaccessory)
unirefaccessorymds = fit(MDS, unirefaccessorydm, distances=true)
unirefaccessorymdsaxes = [v / sum(eigvals(unirefaccessorymds)) for v in eigvals(unirefaccessorymds)]

ukidsunirefaccessorymds = fit(MDS, unirefaccessorydm[ukids,ukids], distances=true)
ukidsunirefaccessorymdsaxes = [v / sum(eigvals(ukidsunirefaccessorymds)) for v in eigvals(ukidsunirefaccessorymds)]

## ## Figure 1A

species_permanovas = runpermanovas(speciesdm, ukids, allmeta, "species")

##
unirefaccessory_permanovas = runpermanovas(unirefaccessorydm, ukids, allmeta, "unirefaccessory")

##
pfamsdm = pairwise(BrayCurtis(), pfams)
pfamsaccessorydm = pairwise(BrayCurtis(), pfamaccessory)
kosdm = pairwise(BrayCurtis(), kos)
kosaccessorydm = pairwise(BrayCurtis(), koaccessory)
ecsdm = pairwise(BrayCurtis(), ecs)
ecsaccessorydm = pairwise(BrayCurtis(), ecaccessory)

pfams_permanovas = runpermanovas(pfamsdm, ukids, allmeta, "pfams")
pfamsaccessory_permanovas = runpermanovas(pfamsaccessorydm, ukids, allmeta, "pfamsaccessory")
kos_permanovas = runpermanovas(kosdm, ukids, allmeta, "kos")
kosaccessory_permanovas = runpermanovas(kosaccessorydm, ukids, allmeta, "kosaccessory")
ecs_permanovas = runpermanovas(ecsdm, ukids, allmeta, "ecs")
ecsaccessory_permanovas = runpermanovas(ecsaccessorydm, ukids, allmeta, "ecsaccessory")


##


allpermanovas = vcat(
    species_permanovas,
    unirefaccessory_permanovas,
    pfams_permanovas,
    pfamsaccessory_permanovas,
    kos_permanovas,
    kosaccessory_permanovas,
    ecs_permanovas,
    ecsaccessory_permanovas
    )
sort!(allpermanovas, [:label, :feature])

r2 = unstack(allpermanovas, :label, :feature, :R2)
select!(r2, ["label", "species", "unirefaccessory", "pfams", "pfamsaccessory", "kos", "kosaccessory", "ecs", "ecsaccessory"])
r2m = Matrix(r2[!,2:end])

q = unstack(allpermanovas, :label, :feature, :q_value)
select!(q, ["label", "species", "unirefaccessory", "pfams", "pfamsaccessory", "kos", "kosaccessory", "ecs", "ecsaccessory"])
qm = Matrix(q[!,2:end])

qa = let M = fill("", size(qm))
    for i in eachindex(qm)
        ismissing(qm[i]) && continue
        if qm[i] < 0.001
            M[i] = "***"
        elseif qm[i] < 0.01
            M[i] = "**"
        elseif qm[i] < 0.1
            M[i] = "*"
        end
    end
    M
end


## Figure 3

abxr = CSV.read("data/uniprot/uniprot-abxr.tsv")
carbs = CSV.read("data/uniprot/uniprot-carbohydrate.tsv")
fa = CSV.read("data/uniprot/uniprot-fa.tsv")
unirefnames = map(u-> match.(r"UniRef90_(\w+)",u).captures[1], featurenames(unirefaccessory))
carbs = Set(carbs.Entry)
carbpos = findall(x-> in(x, carbs), unirefnames)

neuroactive = getneuroactive(unirefnames) # function in accessories.jl

allneuroactive = union([neuroactive[k] for k in keys(neuroactive)]...)

metadatums = [:white_matter_normed,
              :gray_matter_normed,
              :csf_normed,
              :hippocampus_normed,
              :thalamus_normed,
              :corpus_callosum_normed,
              :limbic_normed,
              :subcortex_normed,
              :neocortex_normed,
              :cerebellum_normed]

allfsea = DataFrame(
            geneset   = String[],
            metadatum = String[],
            median    = Float64[],
            pvalue    = Float64[],
            cors      = Vector{Float64}[])

mdcors = Dict(m=>Float64[] for m in metadatums)

for md in metadatums
    @info "Working on $md"
    filt = map(!ismissing, ukidsmeta[!,md])
    cors = cor(ukidsmeta[filt, md], occurrences(view(unirefaccessory, sites=ukids))[:,filt], dims=2)'
    mdcors[md] = filter(!isnan, cors)
    for (key, pos) in pairs(neuroactive)
        allcors = filter(!isnan, cors[pos])
        notcors = filter(!isnan, cors[Not(pos)])
        length(allcors) < 4 && continue
        @info "    $key"
        mwu = MannWhitneyUTest(allcors, notcors)
        m = median(allcors)
        p = pvalue(mwu)
        push!(allfsea, (geneset=key, metadatum=String(md), median=m, pvalue=p, cors=allcors))
    end
end

using StatsBase

ukidsmeta = copy(ukidsmeta)
ukidsmeta[!,:bfnumber] = map(ukidsmeta.breastfeeding) do bf
    ismissing(bf) && return missing
    occursin("formula", bf) && return 0
    occursin("breast", bf) && return 2
    return 1
end
let md = :bfnumber
    @info "Working on $md"
    filt = map(!ismissing, ukidsmeta[!,md])
    cors = Float64[]
    mdcol = disallowmissing(ukidsmeta[filt, md])
    for row in eachrow(occurrences(view(unirefaccessory, sites=ukids))[:,filt])
        row = collect(vec(row))
        push!(cors, corspearman(mdcol, row))
    end

    mdcors[md] = filter(!isnan, cors)
    for (key, pos) in [pairs(neuroactive)...; ["Carbohydrate metabolism"=>carbpos]]
        allcors = filter(!isnan, cors[pos])
        notcors = filter(!isnan, cors[Not(pos)])
        length(allcors) < 4 && continue
        @info "    $key"
        mwu = MannWhitneyUTest(allcors, notcors)
        m = median(allcors)
        p = pvalue(mwu)
        push!(allfsea, (geneset=key, metadatum=String(md), median=m, pvalue=p, cors=allcors))
    end
end

allfsea = DataFrames.transform(groupby(allfsea, :metadatum), :pvalue => (p-> (adjust(collect(p), BenjaminiHochberg()))) => :qvalue)
CSV.write(joinpath(config["output"]["tables"], "allfsea.csv"), allfsea)

## older kids

oldkidsfsea = DataFrame(
            geneset   = String[],
            metadatum = String[],
            median    = Float64[],
            pvalue    = Float64[],
            cors      = Vector{Float64}[])

oldkidsmdcors = Dict(m=>Float64[] for m in metadatums)

for md in metadatums
    @info "Working on $md"
    filt = map(!ismissing, oldkidsmeta[!,md])
    cors = cor(oldkidsmeta[filt, md], occurrences(view(unirefaccessory, sites=oldkidsmeta.sample))[:,filt], dims=2)'
    oldkidsmdcors[md] = filter(!isnan, cors)
    for (key, pos) in pairs(neuroactive)
        oldercors = filter(!isnan, cors[pos])
        notcors = filter(!isnan, cors[Not(pos)])
        length(oldercors) < 4 && continue
        @info "    $key"
        mwu = MannWhitneyUTest(oldercors, notcors)
        m = median(oldercors)
        p = pvalue(mwu)
        push!(oldkidsfsea, (geneset=key, metadatum=String(md), median=m, pvalue=p, cors=oldercors))
    end
end

oldkidsmeta = copy(oldkidsmeta)
oldkidsmeta[!,:bfnumber] = map(oldkidsmeta.breastfeeding) do bf
    ismissing(bf) && return missing
    occursin("formula", bf) && return 0
    occursin("breast", bf) && return 2
    return 1
end
let md = :bfnumber
    @info "Working on $md"
    filt = map(!ismissing, oldkidsmeta[!,md])
    cors = Float64[]
    mdcol = disallowmissing(oldkidsmeta[filt, md])
    for row in eachrow(occurrences(view(unirefaccessory, sites=oldkidsmeta.sample))[:,filt])
        row = collect(vec(row))
        push!(cors, corspearman(mdcol, row))
    end

    oldkidsmdcors[md] = filter(!isnan, cors)
    for (key, pos) in [pairs(neuroactive)...; ["Carbohydrate metabolism"=>carbpos]]
        oldercors = filter(!isnan, cors[pos])
        notcors = filter(!isnan, cors[Not(pos)])
        length(oldercors) < 4 && continue
        @info "    $key"
        mwu = MannWhitneyUTest(oldercors, notcors)
        m = median(oldercors)
        p = pvalue(mwu)
        push!(oldkidsfsea, (geneset=key, metadatum=String(md), median=m, pvalue=p, cors=oldercors))
    end
end

oldkidsfsea = DataFrames.transform(groupby(oldkidsfsea, :metadatum), :pvalue => (p-> (adjust(collect(p), BenjaminiHochberg()))) => :qvalue)

# ## Supplementary Figure 1

function labeldiff(dm, labels)
    u = sort(unique(labels))
    d = Dict(u1 => Dict() for u1 in u)
    for (l1 , l2) in multiset_permutations(repeat(u,2), 2)

        l1pos = findall(isequal(l1), labels)
        l2pos = findall(isequal(l2), labels)
        ds = vec(dm[l1pos, l2pos])
        d[l1][l2] = ds
    end
    d
end


speciesdiffs = labeldiff(speciesdm[ukids, ukids], allmeta.ageLabel[ukids])
unirefaccessorydiffs = labeldiff(unirefaccessorydm[ukids,ukids], allmeta.ageLabel[ukids])
pfamsdiffs = labeldiff(pfamsdm[ukids,ukids], allmeta.ageLabel[ukids])
pfamsaccessorydiffs = labeldiff(pfamsaccessorydm[ukids,ukids], allmeta.ageLabel[ukids])
kosdiffs = labeldiff(kosdm[ukids,ukids], allmeta.ageLabel[ukids])
kosaccessorydiffs = labeldiff(kosaccessorydm[ukids,ukids], allmeta.ageLabel[ukids])
ecsdiffs = labeldiff(ecsdm[ukids,ukids], allmeta.ageLabel[ukids])
ecsaccessorydiffs = labeldiff(ecsaccessorydm[ukids,ukids], allmeta.ageLabel[ukids])



# ## Additional stats

lowidx, upidx = let (l, u) = quantile(skipmissing(ukidsmeta.cogScore), [0.25,0.75])
    lower = findall(s-> !ismissing(s) && s <= l, ukidsmeta.cogScore)
    upper = findall(s-> !ismissing(s) && s >= u, ukidsmeta.cogScore)
    lower, upper
end


quartmeta = copy(ukidsmeta[[lowidx..., upidx...], :])
quartmeta[!,:quartile] = [i <= length(lowidx) ? "bottom 25%" : "top 25%" for i in 1:nrow(quartmeta)]
quartspecies = view(species, sites=quartmeta.sample)
quartspecies = view(quartspecies, species=[sum(row) > 0 for row in eachrow(occurrences(quartspecies))]) |> copy
quartspeciesdm = pairwise(BrayCurtis(), quartspecies)
quartspeciesmds = fit(MDS, quartspeciesdm, distances=true)
quartspeciesmdsaxes = [v / sum(eigvals(quartspeciesmds)) for v in eigvals(quartspeciesmds)]

quartpermanova = permanova(quartspeciesdm, quartmeta.quartile)

MannWhitneyUTest(quartmeta.shannon[1:length(lowidx)], quartmeta.shannon[(length(lowidx)+1):end])
describe(quartmeta.shannon[1:length(lowidx)])
describe(quartmeta.shannon[(length(lowidx)+1):end])

quartiletests = DataFrame()

for (i, row) in enumerate(eachrow(occurrences(quartspecies)))
    row = vec(row)
    sp = featurenames(quartspecies)[i]
    l = row[1:length(lowidx)]
    u = row[(1+length(lowidx)):end]
    mwu = MannWhitneyUTest(l, u)
    assoc = mean(l) > mean(u) ? "-" : "+"
    push!(quartiletests, (
        species=sp,
        association=assoc,
        median_lower = median(l),
        median_upper = median(u),
        nsamples = count(>(0), row),
        pvalue = pvalue(mwu)
    ))
end
filter!(row-> row.nsamples > 30, quartiletests)
quartiletests[!,:qvalue] = adjust(quartiletests.pvalue, BenjaminiHochberg())
sort!(quartiletests, :qvalue)
CSV.write(joinpath(config["output"]["tables"], "quartiletests.csv"), quartiletests)

# ## Older kids

oldkidsspecies = view(species, sites=oldkidsmeta.sample) |> copy

lowidx, upidx = let (l, u) = quantile(skipmissing(oldkidsmeta.cogScore), [0.25,0.75])
    lower = findall(s-> !ismissing(s) && s <= l, oldkidsmeta.cogScore)
    upper = findall(s-> !ismissing(s) && s >= u, oldkidsmeta.cogScore)
    lower, upper
end


quartmeta = copy(oldkidsmeta[[lowidx..., upidx...], :])
quartmeta[!,:quartile] = [i <= length(lowidx) ? "bottom 25%" : "top 25%" for i in 1:nrow(quartmeta)]
quartspecies = view(species, sites=quartmeta.sample)
quartspecies = view(quartspecies, species=[sum(row) > 0 for row in eachrow(occurrences(quartspecies))]) |> copy
quartspeciesdm = pairwise(BrayCurtis(), quartspecies)
quartspeciesmds = fit(MDS, quartspeciesdm, distances=true)
quartspeciesmdsaxes = [v / sum(eigvals(quartspeciesmds)) for v in eigvals(quartspeciesmds)]

quartpermanova = permanova(quartspeciesdm, quartmeta.quartile)
quartiletests = DataFrame()

for (i, row) in enumerate(eachrow(occurrences(quartspecies)))
    row = vec(row)
    sp = featurenames(quartspecies)[i]
    l = row[1:length(lowidx)]
    u = row[(1+length(lowidx)):end]
    assoc = mean(l) > mean(u) ? "-" : "+"
    mwu = MannWhitneyUTest(l, u)
    push!(quartiletests, (
        species=sp,
        association=assoc,
        median_lower = median(l),
        median_upper = median(u),
        nsamples = count(>(0), row),
        pvalue = pvalue(mwu)
    ))
end

filter!(row-> row.nsamples > 30, quartiletests)
quartiletests[!,:qvalue] = adjust(quartiletests.pvalue, BenjaminiHochberg())
sort!(quartiletests, :qvalue)
CSV.write(joinpath(config["output"]["tables"], "oldkidsquartiletests.csv"), quartiletests)

## write transposed file with significant bugs

sigs = filter(row-> row.qvalue < 0.2, quartiletests).species
sigsdf = select(allmeta, [:subject, :timepoint, :sample])
for sp in sigs
    v = vec(occurrences(view(species, sites=sigsdf.sample, species=[sp])))
    sigsdf[:, sp] = v
    allmeta[:, sp] = v
end

CSV.write(joinpath(config["output"]["tables"], "sigcog_species.csv"),sigsdf)



# ## Linear Models
# None of these find anything significant

using GLM

function longtaxfromcomm(cm)
    features = featurenames(cm)
    samples  = samplenames(cm)

    occ = occurrences(cm)

    df = DataFrame((sample=samples[j], taxon=features[i], abundance=occ[i,j])
                        for i in eachindex(features)
                        for j in eachindex(samples))
    return df
end

taxlong = longtaxfromcomm(oldkidsspecies)
taxlong = leftjoin(taxlong, select(oldkidsmeta,
                            [:sample, :childGender, :correctedAgeDays,
                             :mother_HHS, :cogScore, metadatums...]),
                         on=:sample)
taxlong = leftjoin(taxlong, select(quartmeta, :sample, :quartile), on=:sample)
taxlong.asq = asin.(sqrt.(taxlong.abundance))

taxgroups = groupby(taxlong, :taxon)

function runlm!(resultsdf, group, md)
    group = filter(row-> !ismissing(row[md]) && row.abundance > 0, group)
    nrow(group) > 10 || return
    sp = first(group.taxon)
    
    m = @eval lm(@formula(asq ~ $md + correctedAgeDays + childGender + mother_HHS), $group)
    tbl = coeftable(m)
    df = DataFrame([tbl.rownms, tbl.cols...], [:variable, Symbol.(tbl.colnms)...])
    names!(df, [:variable, :estimate, :stderror, :tvalue, :pvalue, :confint5, :confint95])
    df[!, :taxon] .= sp
    append!(resultsdf, df)
end

glms = DataFrame()

for md in metadatums
    @info "GLMs for $md"
    for grp in taxgroups
        runlm!(glms, grp, md)
    end
end

filter!(row-> Symbol(row.variable) in metadatums, glms)

glms = DataFrames.transform(groupby(glms, :variable), :pvalue => (p-> (adjust(collect(p), BenjaminiHochberg()))) => :qvalue)
sort!(glms,:qvalue)

CSV.write(joinpath(config["output"]["tables"], "speciesglms.csv"), glms)

# ## upper/lower quartile

quartglms = DataFrame()

for grp in taxgroups
    grp = filter(row-> !ismissing(row.quartile) && row.abundance > 0, grp)
    nrow(grp) > 10 || continue
    if length(unique(grp.quartile)) == 1
        @info "df for $(first(grp.taxon)) only has 1 quartile ($(first(grp.quartile)))"
        continue
    end

    sp = first(grp.taxon)
    m = lm(@formula(asq ~ quartile + correctedAgeDays + childGender + mother_HHS), grp)
    tbl = coeftable(m)
    df = DataFrame([tbl.rownms, tbl.cols...], [:variable, Symbol.(tbl.colnms)...])
    names!(df, [:variable, :estimate, :stderror, :tvalue, :pvalue, :confint5, :confint95])
    df[!, :taxon] .= sp
    append!(quartglms, df)
end

cogs = findall(row-> startswith(row.variable, "quartile"), eachrow(quartglms))
quartglms.qvalue = Union{Missing,Float64}[missing for _ in 1:nrow(quartglms)]
quartglms.qvalue[cogs] .= adjust(quartglms[cogs,:pvalue], BenjaminiHochberg())
sort!(quartglms, :qvalue)
CSV.write(joinpath(config["output"]["tables"], "quartilespeciesglms.csv"), quartglms)

open(joinpath(config["output"]["other"], "echo_taxa.txt"), "w") do io
    for s in speciesnames(species)
        println(io, s)
    end
end
# ## Exports

using JLD2
@assert sitenames(species) == allmeta.sample
allmeta.pcopri = collect(vec(occurrences(view(species, species=["Prevotella_copri"]))))

@save "analysis/figures/assets/metadata.jld2" allmeta oldkidsmeta ukidsmeta ukids oldkids
@save "analysis/figures/assets/taxa.jld2" species speciesmds speciesmdsaxes ukidsspeciesmds ukidsspeciesmdsaxes
@save "analysis/figures/assets/unirefs.jld2" unirefaccessorymds unirefaccessorymdsaxes ukidsunirefaccessorymds ukidsunirefaccessorymdsaxes
@save "analysis/figures/assets/otherfunctions.jld2" kos kosdiffs kosdm ecs ecsdm ecsdiffs pfams pfamsdiffs pfamsdm koaccessory kosaccessorydiffs kosaccessorydm ecaccessory ecsaccessorydiffs ecsaccessorydm pfamaccessory pfamsaccessorydiffs pfamsaccessorydm
@save "analysis/figures/assets/permanovas.jld2" r2 r2m qa allpermanovas species_permanovas unirefaccessory_permanovas kos_permanovas pfams_permanovas
@save "analysis/figures/assets/fsea.jld2" allfsea oldkidsfsea mdcors oldkidsmdcors
@save "analysis/figures/assets/difs.jld2" speciesdiffs unirefaccessorydiffs kosdiffs pfamsdiffs ecsdiffs  kosaccessorydiffs pfamsaccessorydiffs ecsaccessorydiffs
@save "analysis/figures/assets/stratkos.jld2" stratkos
@save "analysis/figures/assets/cogquartiles.jld2" quartmeta quartspecies quartspeciesdm quartspeciesmds quartspeciesmdsaxes quartiletests

