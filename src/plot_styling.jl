function categorical_colors(items, levels, colors)
    length(levels) <= length(colors) || throw(ArgumentError("Only $(length(colors)) colors for $(length(levels)) levels"))
    colarray = similar(colors, length(items))
    
    colmap = Dict(l => colors[i] for (i,l) in enumerate(levels))
    for i in eachindex(items)
        colarray[i] = colmap[items[i]]
    end
    return colarray
end

function categorical_colors(items, colors)
    levels = unique(items)
    categorical_colors(items, levels, colors)
end

function mds_percent(mds)
    [v / sum(Microbiome.MultivariateStats.eigvals(mds)) for v in Microbiome.MultivariateStats.eigvals(mds)]
end

function mds_percent(mds, inds)
    mds_percent(mds)[inds]
end

mds_format(mds, ind; digits=2) = "MDS $ind ($(round(mds_percent(mds, ind)*100, digits=digits)) %)"