
include("Mesh_utilities.jl")

function remove_outside_tris(p, t, pv)
    pmid = dropdims(sum(p[t,:], dims=2), dims=2) / 3
    is_inside = inpolygon(pmid, pv)
    t = t[is_inside,:]
end

function triarea(p, t)
    d12 = @. p[t[:,2],:] - p[t[:,1],:]
    d13 = @. p[t[:,3],:] - p[t[:,1],:]
    @. abs(d12[:,1] * d13[:,2] - d12[:,2] * d13[:,1]) / 2
end

function remove_tiny_tris(p, t)
    t = t[triarea(p,t) .> 1e-14,:]
end

function circumcenter(p)
    dp1 = @. p[2,:] - p[1,:]
    dp2 = @. p[3,:] - p[1,:]

    mid1 = @. (p[1,:] + p[2,:]) / 2
    mid2 = @. (p[1,:] + p[3,:]) / 2

    s = [ -dp1[2] dp2[2]; dp1[1] -dp2[1]] \ (-mid1 .+ mid2)
    pc = @. mid1' + s[1] * [-dp1[2] dp1[1]]
end

function edge_midpoints(p, t)
    pmid = reshape(p[t,:] + p[t[:,[2,3,1]],:], :, 2) / 2
    pmid = unique(pmid, dims=1)
end

function pmesh(pv, hmax, nref)
    p = zeros(Float64, 0, 2)
    for i = 1:size(pv,1) - 1
        pp = pv[i:i+1,:]
        L = sqrt(sum(diff(pp, dims=1).^2, dims=2))[1]
        if L > hmax
            n = ceil(Int, L / hmax)
            ss = (0:n) / n
            pp = [1 .- ss ss] * pp
        end
        p = [p; pp[1:end-1,:]]
    end

    t = zeros(Int64, 0, 3)
    while true
        t = delaunay(p)
        t = remove_tiny_tris(p, t)
        t = remove_outside_tris(p, t, pv)
        # tplot(p,t), pause(1e-3)

        area = triarea(p, t)
        maxarea, ix = findmax(area)
        if maxarea < hmax^2 / 2
            break
        end
        pc = circumcenter(p[t[ix,:],:])
        p = [p; pc]
    end

    for iref = 1:nref
        p = [p; edge_midpoints(p, t)]
        t = delaunay(p)
        t = remove_tiny_tris(p, t)
        t = remove_outside_tris(p, t, pv)
        # tplot(p, t), pause(1e-3)
    end

    e = boundary_nodes(t)
    p, t, e
end


function fempoi(p,t,e)
    K=size(t)
    N=size(p)
    N=N[1]
    K=K[1]
    c=zeros(3,3)
    Ak=zeros(K,3,3)
    bk=zeros(K)
    for k=1:K
        area=triarea(p,t)
        V=[1 p[t[k,1],1] p[t[k,1],2]; 1 p[t[k,2],1] p[t[k,2],2]; 1 p[t[k,3],1] p[t[k,3],2] ]
        c=V\[1 0 0; 0 1 0; 0 0 1];
        for α=1:3
            for β=1:3
                Ak[k,α,β]=area[k]*(c[2,α]*c[2,β]+c[3,α]*c[3,β])
                bk[k]=area[k]/3;
            end
        end
    end
    A = Tuple{Int64,Int64,Float64}[]  # Array of matrix elements (row,col,value)
    b = zeros(N)

# Main loop, insert stencil in matrix for each node point
    for k=1:K
        for α=1:3
            for β=1:3
                push!(A,(t[k,α],t[k,β],Ak[k,α,β]))
                b[t[k,α]]+=bk[k]
            end
        end
    end
    A = sparse((x->x[2]).(A), (x->x[1]).(A), (x->x[3]).(A), N, N)
    Ne=size(e)
    Ne=Ne[1];
    for i=1:Ne
        for j=1:N
            A[j,e[i,1]]=(e[i,1]==j);
        end
        b[e[i,1]]=0;
    end
    
    return A,b
end

using SparseArrays
using DelimitedFiles
using PyPlot


function gs_iteration(A,b,u)
    p_num = MPI.Comm_size(comm)
    p_id = MPI.Comm_rank(comm)
    u_new = zeros(size(u));
    for row = p_id+1:p_num:size(A,1)
        u_new[row] = b[row];
        aii = 1;
        for j = A.colptr[row]:A.colptr[row+1]-1
            if A.rowval[j]!=row
                u_new[row] = u_new[row] - A.nzval[A.rowval[j]]*u[A.rowval[j]]

            else
                aii = A.nzval[j]
            end
        end
        u_new[row] = u_new[row]/aii;
    end
    return u_new
end

function Gaussian(A,b)
    

    local u = zeros(size(b));
    for t = 1:10000
        u = gs_iteration(A,b,u)
        MPI.Allreduce!(u,+, comm)
    end
    return u 
end
using MPI
MPI.Init()
global comm = MPI.COMM_WORLD
global p_num = MPI.Comm_size(comm)
global p_id = MPI.Comm_rank(comm)
x = 0:.1:1
y = 0.1*(-1).^(0: 10)
pv = [x y; .5 .6; 0 .1]
p, t, e = pmesh(pv, 0.001, 0)
e = e[@. p[e, 2] > (.6 - abs(p[e, 1] - 0.5) - 1e-6) ]

A,b = fempoi(p,t,e)
global N = size(A,1)
start_time = time()
u = Gaussian(A,b)
end_time = time()
if p_id == 0
    writedlm("u.txt",u)
    println((end_time-start_time)/10000)
end


MPI.Barrier(comm)
