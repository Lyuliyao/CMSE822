using PyPlot
using DelimitedFiles
u = readdlm("result.txt",' ')
z = size(u,2);
x = range(0,1,length= z)
plot(x,u')
plot(x,sin.(pi*x))
savefig("result.png")

