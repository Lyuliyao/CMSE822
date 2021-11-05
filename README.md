# CMSE822

In this project, we want to use the MPI and OpenMP to solve following poisson equation. 
$$
\Delta u = f,
$$
subject to appropriate boundary conditions. 

We can use Finite Element method to discrete this problem and solve the system with Gauss seidel iteration, which is easy to parallel computation. 

If time is acceptable, I want to try to solve it on a unstructured mesh (triangle one to be specific) for two dimensional problems, which is more general useful to deal with different structure. There is an potential code in (https://github.com/pgebhardt/libdistmesh)

The main programming languarge will be C++.

The goal of the project is to implement parallel explicit method in 1D and 2D with MPI, OpenMP, and MPI+OpenMP And compare their performance.

The plan for the project is shown blow.

| Time          | Task                                        |
| ------------- | :------------------------------------------ |
| 11.7-11.14    | 1D problem serial version,1D problem OpenMP |
| 11.14 - 11.21 | 1D problem MPI, 1D problem OpenMP+ MPI      |
| 11.21-11.28   | 2D problem serial version,2D problem OpenMP |
| 11.28-12.5    | 2D problem MPI, 2D problem MPI+OpenMP       |
| 12.5-12.12    | finish the report.                          |

