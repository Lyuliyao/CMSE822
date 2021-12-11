using PyPlot
time_ues = [6.73793e-06,4.24195e-06,3.68396e-06,2.9043e-06,3.78522e-06,7.95009e-06]
N_process = [1 2 4 16 32 64];
loglog(N_process',time_ues)
savefig("time_OpenMP_2D.png")