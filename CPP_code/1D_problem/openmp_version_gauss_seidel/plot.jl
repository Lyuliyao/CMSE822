using PyPlot
using DelimitedFiles
time = [0.000306678,0.000166172,7.82919e-05,3.56148e-05,2.05295e-05,1.66777e-05]
thread_num = [1,2,4,8,16,32];
loglog(thread_num,time,LABEL= "number of grid :32000")
time = [0.00141793,0.000603242,0.000272592,0.000132982,4.24693e-05,2.73426e-05]
loglog(thread_num,time,LABEL= "number of grid :64000")
legend()
savefig("open_mp_result.png")

