time_used = [10.25  
5.42
2.96 
1.51 
0.75 
0.38]

thrd_num = [1 
2 
4 
8 
16 
32]
using PyPlot 
figure(figsize=(8,8))
loglog(thrd_num,time_used)

savefig("julia_MPI_2D.png")