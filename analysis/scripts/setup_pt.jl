using PrettyTables

rounder = (v, i, j) -> typeof(v) <: AbstractFloat ? round(v,digits=3) : v
# print ~15 random rows
randrowfilter(data, i) = rand() < (1 / size(data, 1)) * 15
@ptconfclean # clear any previous configuration
@ptconf formatters = rounder nosubheader=true screen_size=(20,120) filters_row=(randrowfilter,)
