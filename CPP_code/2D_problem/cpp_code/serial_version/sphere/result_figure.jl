include("../../../julia_utilities/pmesh.jl")

using PyPlot
using DelimitedFiles
using SparseArrays
u = readdlm("result.txt",' ')
u = u[1:end-1]


n = 32
phi = 2pi*(0: n) /n
pv = [cos.(phi) sin.(phi)]
p, t, e = pmesh(pv, 2pi/n, 0)
fempoi(p, t, e)
tplot(p, t, u)

savefig("result.png")