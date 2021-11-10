#in this code we test the 1D serial version of the finite element method
# the basis function hat function is used.

function mass_matric(N::Integer)
    # here N denote divide the line into N part.
    M = zeros(N+1,N+1);
    for i = 2:N
        M[i,i] = 2*N;
        M[i,i-1] = -N;
        M[i,i+1] = -N;
    end
    return M
end

function fcn_left(x)
    return pi^2*sin(pi*x)
end

function poisson_1D(Dbc_value)
    N = 10;
    M = mass_matric(N)
    M[1,1] = 1;

    M[N+1,N+1] = 1;
    L = zeros(N+1,1);
    L[1] = Dbc_value[1]
    L[end] = Dbc_value[2]
    f_l(x,n) = (x+1)/N/4*fcn_left((x+1)/2/N+(n-1)/N)
    f_r(x,n) = (1-(x+1)/2)/N/2*fcn_left((x+1)/N/2+n/N)
    for i = 2:N
        L[i] = f_l(-1/sqrt(3),i-1) + f_l(1/sqrt(3),i-1)
        L[i] = L[i] + f_r(-1/sqrt(3),i-1) + f_r(1/sqrt(3),i-1)
    end
    u[1] = Dbc_value[1]
    L[end] = Dbc_value[2]
    for i = 2:N
        u[i] = 1/M[i,i]*(L[i] - M[i,i-1]*u[i-1] - M[i,i+1] *u[i+1])
    end
    u = M\L
end
u = poisson_1D([0,0])
plot(u)
