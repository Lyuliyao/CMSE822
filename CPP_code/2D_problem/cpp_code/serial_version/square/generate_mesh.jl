using SparseArrays
using DelimitedFiles
using PyPlot
include("../../../julia_utilities/pmesh.jl")
pv = Float64[0 0; 1 0; 1 1; 0 1; 0 0]
p,t,e = pmesh(pv, 0.15, 0)
e = e[@. (p[e, 1] < 1e-6) | (p[e, 2] < 1e-6) ]
fempoi(p, t, e)


