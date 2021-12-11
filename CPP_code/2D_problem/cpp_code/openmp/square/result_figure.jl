include("../../../julia_utilities/pmesh.jl")

using PyPlot
using DelimitedFiles
using SparseArrays
u = readdlm("result.txt",' ')
u = u[1:end-1]


pv = Float64[0 0; 1 0; 1 1; 0 1; 0 0]
p,t,e = pmesh(pv, 0.15, 0)
e = e[@. (p[e, 1] < 1e-6) | (p[e, 2] < 1e-6) ]
tplot(p, t, u)

savefig("result.png")