# one dimensional Legendre polynomials

@inline function getLegendrePn(x, n::Int64)
    #= Compute Pn evaluated at x ∈ [-1,1] =#
    if n>1
        resm = ones(Float64, length(x))
        res  = copy(x)
        resp = zeros(Float64, length(x))
        for i=2:n
            resp .= ((2*i - 1) * x .* res - (i-1)*resm) / i
            resm .= res; res .= resp
        end
        return res
    elseif n==0
        return ones(Float64, length(x))
    elseif n==1
        return x
    end
end

@inline function getAllLegendrePn(x, n::Int64)
    #= Compute all Pn for all x ∈ [-1, 1] =#
    nx = length(x)
    if n > 1
        res = zeros(Float64, nx, n+1)
        res[:, 2] .= x
        for i=2:n
            res[:, i+1] = ( (2*i -1) * x .* res[:, i] - (i-1)*res[:, i-1]) / i
        end
        return res
    elseif n == 0
        return zeros(Float64, nx, 1)
    elseif n == 1
        res = zeros(Float64, nx, 2)
        res[:, 2] .= x
        return res
    end
end

# interior monomials
function getMonomialsP1(xE::Array{Float64, 1}, hE::Float64)
    #=
    Return an array of Λ-functions enabling the computation of the
    P1 monomials
    =#
    return [x -> 1., x -> (x[1] - xE[1])/hE, x -> (x[2] - xE[2])/hE]
end

function getMonomialsP1(xE::Array{Float64, 1}, hE::Float64, Cs::Array{Float64,1})
    #= Compute single P1 monomials with 0 average =#
    return [x -> 1., x -> (x[1] - xE[1])/hE + Cs[2], x -> (x[2] - xE[2])/hE + Cs[3]]
end

function getGradientMonomialsP1(xE::Array{Float64, 1}, hE::Float64)
    #= Compute the gradient of the single P1 monomials =#
    return [x -> [0.,0.], x -> [1/hE, 0.], x -> [0., 1/hE]]
end

function getMonomialsP2(xE::Array{Float64, 1}, hE::Float64)
    #= Compute single P2 monomials =#
    hE2 = hE^2
    return  [
        x -> 1,
        x -> (x[1] - xE[1])/hE,
        x -> (x[2] - xE[2])/hE,
        x -> (x[1] - xE[1])^2 / hE2,
        x -> (x[1] - xE[1])*(x[2] - xE[2])/hE2,
        x -> (x[2] - xE[2])^2 / hE2]
end

function getMonomialsP2(xE::Array{Float64, 1}, hE::Float64, Cs::Array{Float64,1})
    #= Compute single P2 monomials with 0 average =#
    hE2 = hE^2
    return [
        x -> 1,
        x -> (x[1] - xE[1])/hE + Cs[2],
        x -> (x[2] - xE[2])/hE + Cs[3],
        x -> (x[1] - xE[1])^2 / hE2 + Cs[4],
        x -> (x[1] - xE[1])*(x[2] - xE[2])/hE2 + Cs[5],
        x -> (x[2] - xE[2])^2 / hE2 + Cs[6]
    ]
end

function getGradientMonomialsP2(xE::Array{Float64, 1}, hE::Float64)
    #= Compute gradient of the single P2 monomials =#
    hE2 = hE^2
    return [
        x -> [0.,0.],
        x -> [1/hE, 0.],
        x -> [0., 1/hE],
        x -> [2*(x[1] - xE[1])/hE2, 0.],
        x -> [(x[2] - xE[2])/hE2, (x[1] - xE[1])/hE2],
        x -> [0., 2*(x[2] - xE[2])/hE2]
    ]
end
