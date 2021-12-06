using PyPlot
time = [0.000365184,0.000250849,0.000217986,0.00023789,0.000324083,0.000285311]
thread_num = [1,2,4,8,16,32];
for i in time
      @printf("%1.2e &",i)
end
loglog(thread_num,time,label= "number of grid :32000")
time = [0.000804676,0.000510317,0.000434813,0.000435639,0.000501324,0.000701634]
for i in time
      @printf("%1.2e &",i)
end
loglog(thread_num,time,label= "number of grid :64000")
legend()
xlabel("numer of thread")
ylabel("time")
savefig("mpi_result.png")
