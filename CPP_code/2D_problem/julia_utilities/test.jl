
using SparseArrays
using DelimitedFiles
using PyPlot
include("pmesh.jl")
figure(1)
pv = Float64[0 0; 1 0; 1 1; 0 1; 0 0]
p,t,e = pmesh(pv, 0.15, 0)
e = e[@. (p[e, 1] < 1e-6) | (p[e, 2] < 1e-6) ]
fempoi(p, t, e)
tplot(p, t)
savefig("square.png")
figure(2)
x = 0:pi/5:2*pi
pv = zeros(Float64,size(x,1),2)
pv[:,1] = sin.(x)
pv[:,2] = cos.(x)
p,t,e = pmesh(pv, 0.15, 0)
e = e[@. (p[e, 1] < 1e-6) | (p[e, 2] < 1e-6) ]
fempoi(p, t, e)
tplot(p, t)
savefig("five_line.png")
figure(3)
x = 0:pi/24:2*pi
pv = zeros(Float64,size(x,1),2)
pv[:,1] = sin.(x)
pv[:,2] = cos.(x)
p,t,e = pmesh(pv, 0.15, 0)
e = e[@. (p[e, 1] < 1e-6) | (p[e, 2] < 1e-6) ]
fempoi(p, t, e)
tplot(p, t)
savefig("sphere.png")
figure(4)
x = 0:0.1:2*pi/3
pv = zeros(Float64,size(x,1),2)
pv[:,1] = sin.(x)
pv[:,2] = cos.(x)
pv[1,1] = 0
pv[1,2] = 0
pv[end,1] = 0
pv[end,2] = 0
p,t,e = pmesh(pv, 0.15, 0)
e = e[@. (p[e, 1] < 1e-6) | (p[e, 2] < 1e-6) ]
fempoi(p, t, e)
p = tplot(p, t)
savefig("half_sphere.png")
