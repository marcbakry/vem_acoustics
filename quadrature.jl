#=
This file contains all of the functions to perform the integration on triangles.
=#

mutable struct Quadrature1d
    wg::Array{Float64, 1}
    xg::Array{Float64, 1}
end

mutable struct Quadrature2d
    wg::Array{Float64, 1}
    xg::Array{Float64, 2}
end

Quadrature1d() = Quadrature1d(zeros(Float64, 0), zeros(Float64, 0))
Quadrature2d() = Quadrature2d(zeros(Float64, 0), zeros(Float64, 0, 0))

function quadGL1d(ng::Int64)
    #= Hard-coded quadrature rules on [-1,1] =#
    if ng == 1
        xg = [0.0]; wg = [2.0];
    elseif ng == 2
        xg = [-0.577350269189625764509148780502, 0.577350269189625764509148780502]; wg = [1.0, 1.0]
    elseif ng == 3
        xg = [-0.774596669241483377035853079956, 0.0, 0.774596669241483377035853079956]; wg = [5. /9., 8. /9., 5. /9.]
    elseif ng == 4
        xg = [-0.861136311594052575223946488893, -0.339981043584856264802665759103, 0.339981043584856264802665759103, 0.861136311594052575223946488893]; wg = [0.347854845137453857373063949222, 0.652145154862546142626936050778, 0.652145154862546142626936050778,0.347854845137453857373063949222];
    elseif ng == 5
        xg = [-0.906179845938663992797626878299, -0.538469310105683091036314420700, 0.0, 0.538469310105683091036314420700, 0.906179845938663992797626878299]; wg = [0.236926885056189087514264040720, 0.478628670499366468041291514836, 0.568888888888888888888888888889, 0.478628670499366468041291514836, 0.236926885056189087514264040720];
    else
        xg = []; wg = [];
    end
    return Quadrature1d(wg, xg)
end

function quadGL1d(ng::Int64, a::Float64, b::Float64)
    #= Quadrature rule on [a,b] =#
    ref = quadGL1d(ng)
    cw = (b - a)/2.
    cx = (a + b)/2.
    for i = 1:ng
        ref.wg[i] *= cw
        ref.xg[i] =  cx + cw*ref.xg[i]
    end
    return ref
end

function quadGL1d(ref::Quadrature1d, a::Float64, b::Float64)
    #= Same as before, but from an already-existing rule =#
    cw = (b - a)/2.
    cx = (a + b)/2.
    xg = copy(ref.xg)
    wg = copy(ref.wg)
    for i = 1:length(ref.wg)
        wg[i] *= cw
        xg[i] =  cx + cw*xg[i]
    end
    return Quadrature1d(wg, xg)
end

function quadGL2d(ng::Int64)
    #= Hard-coded 2d quadrature rules on {(0,0),(1,0),(0,1)} =#
    if ng == 1
        xg = reshape([1/3.; 1/3.], (2, 1)); wg = [0.5]
    elseif ng == 3
        xg = [2/3 1/6 1/6; 1/6 2/3 1/6]; wg = [1/6., 1/6., 1/6.];
    elseif ng == 4
        xg = [1/3 0.6 0.2 0.2; 1/3 0.2 0.6 0.2]; wg = [-27/96., 25. /96., 25. /96., 25. /96.];
    elseif ng == 6
        a1 = 0.816847572980459; b1 = 0.091576213509771; w1 = 0.5*0.109951743655322;
        a2 = 0.108103018168070; b2 = 0.445948490915965; w2 = 0.5*0.223381589678011;
        xg = [a1 b1 b1 a2 b2 b2; b1 a1 b1 b2 a2 b2]; wg = [w1, w1, w1, w2, w2, w2];
    else
        xg = []; wg = [];
    end
    return Quadrature2d(wg, xg)
end

function quadGL2d(ng::Int64, T::Array{Float64, 2})
    #= Returns quadrature on T =#
    ref = quadGL2d(ng)
    A  = T[:, 1]
    v1 = T[:, 2] - T[:, 1]
    v2 = T[:, 3] - T[:, 1]
    ac = abs(v1[1]*v2[2] - v1[2]*v2[1])
    xg = zeros(Float64, 2, ng)
    wg = zeros(Float64, ng)
    for i=1:ng
        xg[1, i] = A[1] + v1[1]*ref.xg[1, i] + v2[1]*ref.xg[2, i]
        xg[2, i] = A[2] + v1[2]*ref.xg[1, i] + v2[2]*ref.xg[2, i]
        wg[i]    = ac * ref.wg[i]
    end
    return Quadrature2d(wg, xg)
end

function quadGL2d(ref::Quadrature2d, T::Array{Float64, 2})
    #= Quadrature on triangle T from a reference quadrature ref =#
    ng = length(ref.wg)
    A  = T[:, 1]
    v1 = T[:, 2] - T[:, 1]
    v2 = T[:, 3] - T[:, 1]
    ac = abs(v1[1]*v2[2] - v1[2]*v2[1])
    xg = zeros(Float64, 2, ng)
    wg = zeros(Float64, ng)
    for i=1:ng
        xg[1, i] = A[1] + v1[1]*ref.xg[1, i] + v2[1]*ref.xg[2, i]
        xg[2, i] = A[2] + v1[2]*ref.xg[1, i] + v2[2]*ref.xg[2, i]
        wg[i]    = ac * ref.wg[i]
    end
    return Quadrature2d(wg, xg)
end

function integrate(quad::Quadrature1d, f::Function)
    #= Computes the integral of f using quad =#
    ng  = length(quad.wg)
    tot = zero(f(0.))
    for n=1:ng
        @inbounds tot += quad.wg[n] * f(quad.xg[n])
    end
    return tot
end

function integrate(quad::Quadrature2d, f::Function)
    #= same as before, but in 2d =#
    ng  = length(quad.wg)
    tot = zero(f([0., 0.]))
    for n=1:ng
        @inbounds tot += quad.wg[n] * f(quad.xg[:, n])
    end
    return tot
end
