include("../../../julia_utilities/pmesh.jl")

using PyPlot
using DelimitedFiles
using SparseArrays
u = readdlm("result.txt",' ')
u = u[1:end-1]


p = read("p.txt")

t = read("t.txt")

tplot(p, t, u)

savefig("result.png")