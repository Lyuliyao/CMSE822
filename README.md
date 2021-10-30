# CMSE822

In this project, we want to use the MPI and OpenMP to solve following poisson equation. 
$$
\Delta u = f,
$$
subject to appropriate boundary conditions. 

We can use Finite Element method to discrete this problem and solve the system with Gauss seidel iteration, which is easy to parallel computation. 

If time is acceptable, I want to try to solve it on a unstructured mesh (triangle one to be specific) which is more general useful to deal with different structure. There is an potential code in (https://github.com/pgebhardt/libdistmesh)

I can also implement it in Julia (https://docs.julialang.org/en/v1/manual/distributed-computing/), which I believe is an promising language, even though I'm not sure whether the multiprocessing part of it have been well developed or not.


The goal of the project is to implement parallel explicit method in 1D and 2D with MPI, OpenMP, and MPI+OpenMP. And compare their performance.