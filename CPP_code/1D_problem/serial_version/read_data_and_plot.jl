using PyPlot
using DelimitedFiles
u = readdlm("result.txt",' ')

plot([0:1/11:1],u)
