
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
        for ??=1:3
            for ??=1:3
                Ak[k,??,??]=area[k]*(c[2,??]*c[2,??]+c[3,??]*c[3,??])
                bk[k]=area[k]/3;
            end
        end
    end
    A = Tuple{Int64,Int64,Float64}[]  # Array of matrix elements (row,col,value)
    b = zeros(N)

# Main loop, insert stencil in matrix for each node point
    for k=1:K
        for ??=1:3
            for ??=1:3
                push!(A,(t[k,??],t[k,??],Ak[k,??,??]))
                b[t[k,??]]+=bk[k]
            end
        end
    end
    A = sparse((x->x[1]).(A), (x->x[2]).(A), (x->x[3]).(A), N, N)
    Ne=size(e)
    Ne=Ne[1];
    for i=1:Ne
        for j=1:N
            A[e[i,1],j]=(e[i,1]==j);
        end
        b[e[i,1]]=0;
    end
    u = A\b
    return u
end

using SparseArrays
using DelimitedFiles
using PyPlot

n = 32
phi = 2pi*(0: n) /n
pv = [cos.(phi) sin.(phi)]
p, t, e = pmesh(pv, 2pi/n, 0)
u = fempoi(p, t, e)
tplot(p, t, u)