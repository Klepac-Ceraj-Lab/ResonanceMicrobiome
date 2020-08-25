const color1 = ColorSchemes.Set1_9.colors
const color2 = ColorSchemes.Set2_8.colors
const color3 = ColorSchemes.Set3_12.colors
const color4 = ColorSchemes.Paired_12.colors

function get_neuroactive_kos(neuroactivepath="data/uniprot/gbm.txt")
    neuroactive = HashDictionary{String, Vector{String}}()
    desc = ""
    for line in eachline(neuroactivepath)
       line = split(line, r"[\t,]")
       if startswith(line[1], "MGB")
           (mgb, desc) = line
           desc = rstrip(replace(desc, r"\bI+\b.*$"=>""))
           desc = replace(desc, r" \([\w\s]+\)$"=>"")
           desc = replace(desc, r"^.+ \(([\w\-]+)\) (.+)$"=>s"\1 \2")
           @info "getting unirefs for $desc"
           !in(desc, keys(neuroactive)) && insert!(neuroactive, desc, String[])
       else
           filter!(l-> occursin(r"^K\d+$", l), line)
           append!(neuroactive[desc], String.(line))
       end
   end
   return neuroactive
end


function getneuroactive(features, neuroactivepath="data/uniprot/gbm.txt", mappath="/babbage/echo/profiles/map_ko_uniref90.txt.gz")
    neuroactivekos = get_neuroactive_kos(neuroactivepath)

    kos2uniref = Dictionary{String, Vector{String}}()
    open(mappath) do io
        for line in eachline(GzipDecompressorStream(io))
            line = split(line, '\t')
            insert!(kos2uniref, line[1], map(x-> String(match(r"UniRef90_(\w+)", x).captures[1]), line[2:end]))
        end
    end

    neuroactive_index = Dictionary{String, Vector{Int}}()
    for na in keys(neuroactivekos)
        searchfor = Iterators.flatten([kos2uniref[ko] for ko in neuroactivekos[na] if ko in keys(kos2uniref)]) |> Set
        pos = findall(u-> u in searchfor, features)
        insert!(neuroactive_index, na, pos)
    end
    for k in keys(neuroactive_index)
        unique!(neuroactive_index[k])
    end
    return neuroactive_index
end

fsea(cors, pos) = (cors, pos, MannWhitneyUTest(cors[pos], cors[Not(pos)]))

function fsea(cors, allfeatures::AbstractVector, searchset::Set)
    pos = findall(x-> x in searchset, allfeatures)

    return fsea(cors, pos)
end

function fsea(occ::AbstractMatrix, metadatum::AbstractVector, pos::AbstractVector{<:Int})
    let notmissing = map(!ismissing, metadatum)
        occ = occ[:, notmissing]
        metadatum = metadatum[notmissing]
    end

    cors = cor(metadatum, occ, dims=2)'
    return fsea(cors, pos)
end


function accessorygenes(cm, calcs; lower=0., upper=0.95)
    prev, acc = Bool[], Bool[]

     @showprogress map(eachrow(occurrences(cm))) do row
        c = [calc(row) for calc in calcs]
        push!(prev, any(>(lower), c))
        push!(acc,  all(<(upper), c))
    end
    return view(cm, species=prev), view(cm, species=acc)
end

function runpermanovas(dm, ufilter, md)
    udm = dm[ufilter, ufilter]
    umeta = view(md, ufilter, :)
    perms = vcat(
        permanova(dm, [ismissing(x) ? missing : string(x) for x in md.subject], label="subject"),
        permanova(udm, udm.correctedAgeDays, label="age"),
        permanova(udm, udm,
            datafilter=row-> !ismissing(row.ageLabel) && row.ageLabel != "1 and under",
            fields=[:correctedAgeDays], label="1+ age"),
        permanova(udm, [ismissing(x) ? missing : string(x) for x in udm.birthType], label="birth type"),
        permanova(udm, [ismissing(x) ? missing : string(x) for x in udm.childGender], datafilter=x-> x != "Don't know", label="gender"),
        permanova(udm, udm.mother_HHS, label="mother SES"),
        permanova(udm, udm, fields=[:correctedAgeDays,:limbic_normed], label="limbic volume")[2:2,:],
        permanova(udm, udm, fields=[:correctedAgeDays,:subcortical_normed], label="subcortical volume")[2:2,:],
        permanova(udm, udm, fields=[:correctedAgeDays,:neocortical_normed], label="neocortical volume")[2:2,:],
        permanova(udm, udm, fields=[:correctedAgeDays,:cerebellar_normed], label="cerebellar volume")[2:2,:],
        permanova(udm, udm.cogScore, label="cognitive function"),
        permanova(udm, [ismissing(x) ? missing : string(x) for x in udm.breastfeeding], label="breastfeeding"),
        permanova(udm, [ismissing(x) ? missing : string(x) for x in udm.simple_race], label="race"),
        permanova(udm, udm.childBMI, label="BMI")
        )

    filter!(r-> !ismissing(r[Symbol("Pr(>F)")]), perms)
    perms[!, :feature] .= "species"
    rename!(perms, Symbol("Pr(>F)")=>:p_value)
    disallowmissing!(perms)
    perms.q_value = adjust(perms.p_value, BenjaminiHochberg())
    sort!(perms, :q_value)
end