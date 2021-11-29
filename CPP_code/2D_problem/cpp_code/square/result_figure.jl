include("../../julia_utilities/pmesh.jl")

using PyPlot
using DelimitedFiles
u = readdlm("result.txt",' ')

tplot(p, t, u)

