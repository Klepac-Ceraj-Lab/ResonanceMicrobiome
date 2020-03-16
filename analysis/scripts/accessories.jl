using ColorBrewer
using HypothesisTests
using InvertedIndices
using Statistics
using Dictionaries

const color1 = ColorBrewer.palette("Set1", 9)
const color2 = ColorBrewer.palette("Set2", 8)
const color3 = ColorBrewer.palette("Set3", 12)
const color4 = ColorBrewer.palette("Paired", 12)

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

function getneuroactive(features, neuroactivepath="data/uniprot/gbm.txt")
    neuroactivekos = get_neuroactive_kos(neuroactivepath)

    kos2uniref = Dict()
    for line in eachline("data/engaging/ko2uniref90.txt")
        line = split(line, '\t')
        kos2uniref[line[1]] = map(x-> String(match(r"UniRef90_(\w+)", x).captures[1]), line[2:end])
    end

    neuroactive_index = HashDictionary{String, Vector{Int}}()
    for na in keys(neuroactivekos)
        searchfor = [kos2uniref[ko] for ko in neuroactivekos[na] if ko in keys(kos2uniref)]
        pos = findall(u-> u in searchfor, features)
        insert!(neuroactive, na, pos)
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
