using ColorBrewer
using HypothesisTests
using InvertedIndices
using Statistics
using ProgressMeter

# const color1 = ColorBrewer.palette("Set1", 9)
# const color2 = ColorBrewer.palette("Set2", 8)
# const color3 = ColorBrewer.palette("Set3", 12)
# const color4 = ColorBrewer.palette("Paired", 12)

function getneuroactive(features, neuroactivepath="data/uniprot/gbm.txt")

    kos2uniref = Dict()
    for line in eachline("data/engaging/ko2uniref90.txt")
        line = split(line, '\t')
        kos2uniref[line[1]] = map(x-> String(match(r"UniRef90_(\w+)", x).captures[1]), line[2:end])
    end


    neuroactive = Dict()
    desc=""
    for line in eachline("data/uniprot/gbm.txt")
        line = split(line, r"[\t,]")
        if startswith(line[1], "MGB")
            (mgb, desc) = line
            desc = rstrip(replace(desc, r"\bI+\b.*$"=>""))
            desc = replace(desc, r" \([\w\s]+\)$"=>"")
            desc = replace(desc, r"^.+ \(([\w\-]+)\) (.+)$"=>s"\1 \2")
            @info "getting unirefs for $desc"
            if !in(desc, keys(neuroactive))
                neuroactive[desc] = Int[]
            end
        else
            filter!(l-> occursin(r"^K\d+$", l), line)
            searchfor = [kos2uniref[k] for k in line if k in keys(kos2uniref)]
            length(searchfor) == 0 && continue
            searchfor = Set(vcat(searchfor...))
            pos = findall(u-> u in searchfor, features)
            append!(neuroactive[desc], pos)
        end
    end
    for k in keys(neuroactive)
        unique!(neuroactive[k])
    end
    return neuroactive
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

#
# function plotfsea(df::AbstractDataFrame, ymin, ymax)
#     md = first(eachrow(df)).metadatum
#     plt = boxplot(df.cors, color=:lightgrey,
#         legend=false, xticks=(1:nrow(df), df.geneset), xrotation=90,
#         ylabel="Median $md Correlation",
#         yticks= ymin:0.1:ymax)
#     for (i, row) in enumerate(eachrow(df))
#         if row.qvalue < 0.001
#             annotate!(i, ymax -0.1, "***", align=:center)
#         elseif row.qvalue < 0.01
#             annotate!(i, ymax -0.1, "**", align=:center)
#         elseif row.qvalue < 0.1
#             annotate!(i, ymax -0.1, "*", align=:center)
#         end
#     end
#     return(plt)
# end
