using SparseArrays
using DelimitedFiles
using PyPlot
include("../../../julia_utilities/pmesh.jl")
n = 32
phi = 2pi*(0: n) /n
pv = [cos.(phi) sin.(phi)]
p, t, e = pmesh(pv, 2pi/n, 0)
fempoi(p, t, e)


