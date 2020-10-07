# utilitaries
using LinearAlgebra

include("quadrature.jl")
# include("connectique.jl")

@inline function mod_wrap(x, a)
    #= return the index x modulo a =#
    return mod(x-1, a) + 1
end

function polygonArea(vtx)
    area = 0.
    N = size(vtx, 2)
    for i=1:N
        j = mod_wrap(i+1, N)
        area += vtx[1, i]*vtx[2, j] - vtx[2, i]*vtx[1, j]
    end
    return 0.5*abs(area)
end

function polygonDiameter(vtx)
    diameter = 0.
    N = size(vtx, 2)
    for i=1:(N-1)
        for j=(i+1):N
            diameter = max(diameter, norm(vtx[:, i] - vtx[:, j]))
        end
    end
    return diameter
end

function polygonCentroid(vtx)
    N        = size(vtx, 2)
    centroid = zeros(Float64, 2)
    #
    den = 0.
    for i=1:N
        j = mod_wrap(i+1, N)
        w = (vtx[1, i]*vtx[2, j] - vtx[1, j]*vtx[2, i])
        centroid[1] += (vtx[1, i] + vtx[1, j])*w
        centroid[2] += (vtx[2, i] + vtx[2, j])*w
        den += w
    end
    #
    return centroid/(3 * den)
end

# INTEGRATION OF MULTIPLE FUNCTIONS OVER A POLYGON
function integratePolygon(pverts::Array{Float64,2}, xE::Array{Float64,1},
    ref::Quadrature2d, fs::Array{Function, 1})
    #=
    Integration of all functions in fs over a polygon by subdividing it
    into triangles
    =#
    n   = size(pverts, 2)
    nf  = length(fs)
    res = zeros(Float64, nf)
    for i=1:n
        ip = mod_wrap(i+1, n)
        @inbounds quadT = quadGL2d(ref, [pverts[1,i] pverts[1,ip] xE[1];
        pverts[2,i] pverts[2,ip] xE[2]])
        for f = 1:nf
            @inbounds res[f] += integrate(quadT, fs[f])
        end
    end
    return  res
end


function integratePolygon(pverts::Array{Float64,2}, xE::Array{Float64,1},
    ref::Quadrature2d, fs::Function)
    #= Same as before, but for a single function fs =#
    n   = size(pverts, 2)
    res = 0.
    for i=1:n
        ip = mod_wrap(i+1, n)
        @inbounds quadT = quadGL2d(ref, [pverts[1,i] pverts[1,ip] xE[1];
        pverts[2,i] pverts[2,ip] xE[2]])
        res += integrate(quadT, fs)
    end
    return  res
end
