include("../../../julia_utilities/pmesh.jl")

using PyPlot
using DelimitedFiles
using SparseArrays
u = readdlm("result.txt",' ')
u = u[1:end-1]


x = 0:.1:1
y = 0.1*(-1).^(0: 10)
pv = [x y; .5 .6; 0 .1]
p, t, e = pmesh(pv, 0.04, 0)
e = e[@. p[e, 2] > (.6 - abs(p[e, 1] - 0.5) - 1e-6) ]
tplot(p, t, u)

savefig("result.png")