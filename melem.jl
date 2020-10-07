include("fbase.jl")
include("quadrature.jl")
include("utils.jl")

include("melem_rt0.jl")

using LinearAlgebra

function melemk0P2(verts::Array{Float64, 2})
    #= Compute the A = ∫div(v) div(u) and the B = ∫v u element-matrices when the
        approximation order is k=0. In that case, there are no interior degrees
        of freedom and the amount of dofs is exactly the number of edges =#
    nvtx = size(verts, 2)
    ndof = nvtx
    # 0  : preliminaries
    E  = polygonArea(verts)         # area
    xE = polygonCentroid(verts)     # centroid
    hE = polygonDiameter(verts)     # diameter

    # 1st: compute A; it is easily done thanks to k=0
    DIV = ones(Float64, ndof, ndof) / E

    # 2nd: compute B; somehow a little more difficult as a projection needs to
    # be computed
    gmα = getGradientMonomialsP2(xE, hE)
    ref2d = quadGL2d(3)

    Παβ = zeros(Float64, 5, 5)
    for α=1:5
        for β = 1:5
            Παβ[α,β] = integratePolygon(verts, xE, ref2d,
                x -> dot(gmα[α+1](x), gmα[β+1](x)))
        end
    end
    # RHS: first term: null!
    # RHS: second term
    # initialization of the interior basis functions
    mα = getMonomialsP2(xE, hE)
    Cs = (-1/E)*integratePolygon(verts, xE, ref2d, mα) # average over polygon
    mα = getMonomialsP2(xE, hE, Cs) # 0-averaged monomials
    #
    rhs = zeros(Float64, 5, ndof)
    ref1d = quadGL1d(2, 0., 1.)
    # rhs[α, i] = int_{e_i}{((φ_i⋅n) Φ_α)} where φ_i⋅n = 1
    for i=1:ndof
        ip  = mod_wrap(i+1, ndof)        # return the next index modulo ndof
        e_v = verts[:, ip] - verts[:, i] # tangent to edge
        le  = norm(e_v)                  # length of edge
        for α=1:5
            # integration with quadrature rule
            rhs[α, i] = integrate(ref1d, t -> mα[α+1](verts[:, i] + t*e_v)/le)*le
        end
    end
    Π = Παβ \ rhs # Π is the projection matrix
    # now the Di,α = dof_i(m_α) = ∫(∇mα . n) q^i  matrix
    D = zeros(Float64, ndof, 5)
    for i=1:ndof
       ip  = mod_wrap(i+1, ndof)
       e_v = verts[:, ip] - verts[:, i]
       le  = norm(e_v)
       n_v = e_v / le          # normalized tangent
       n_v = [n_v[2], -n_v[1]] # outside pointing normal vector
       for α = 1:5
           D[i,α] = integrate(ref1d,
            t -> dot(gmα[α+1]( verts[:, i] + t*e_v ), n_v))*le
       end
    end
    # stability constant
    σ    = 1.
    F    = Matrix(1.0I, ndof, ndof) - D*Π
    MASS = Π'*(Παβ*Π) + σ*(F' * F)

    # 3rd:  both parts have been computed
    return DIV, MASS
end


function melemk0P1(verts::Array{Float64, 2})
    #= Compute the A = ∫div(v) div(u) and the B = ∫v u element-matrices when the
        approximation order is k=0. In that case, there are no interior degrees
        of freedom and the amount of dofs is exactly the number of edges =#
    nvtx = size(verts, 2)
    ndof = nvtx
    # 0  : preliminaries
    E  = polygonArea(verts)         # area
    xE = polygonCentroid(verts)     # centroid
    hE = polygonDiameter(verts)     # diameter

    # 1st: compute A; it is easily done thanks to k=0
    DIV = ones(Float64, ndof, ndof) / E

    # 2nd: compute B; somehow a little more difficult as a projection needs to
    # be computed
    gmα = getGradientMonomialsP1(xE, hE)
    ref2d = quadGL2d(1)

    Παβ = zeros(Float64, 2, 2)
    for α=1:2
        for β = 1:2
            # quadrature rule not really necessary as everythin is constant
            # but it mimics the code of melemk0P2
            Παβ[α,β] = integratePolygon(verts, xE, ref2d,
                x -> dot(gmα[α+1](x), gmα[β+1](x)))
        end
    end
    # RHS: first term: null!
    # RHS: second term
    # initialization of the interior basis functions
    mα = getMonomialsP1(xE, hE)
    Cs = (-1/E)*integratePolygon(verts, xE, ref2d, mα) # average over polygon
    mα = getMonomialsP1(xE, hE, Cs) # 0-averaged monomials
    #
    rhs = zeros(Float64, 2, ndof)
    ref1d = quadGL1d(2, 0., 1.)
    # rhs[α, i] = int_{e_i}{((φ_i⋅n) Φ_α)} where φ_i⋅n = 1/le
    for i=1:ndof
        ip  = mod_wrap(i+1, ndof)
        e_v = verts[:, ip] - verts[:, i]
        le  = norm(e_v)
        for α=1:2
            # rhs[α, i] = integrate(ref1d, t -> mα[α+1](verts[:, i] + t*e_v)/le)*le
            rhs[α, i] = mα[α+1]((verts[:, ip] + verts[:, i])/2); # with the midpoint rule
        end
    end
    Π = Παβ \ rhs # Π is the projection matrix
    # now the Di,α = dof_i(m_α) = ∫(∇mα . n) q^i  matrix
    D = zeros(Float64, ndof, 2)
    for i=1:ndof
       ip  = mod_wrap(i+1, ndof)
       e_v = verts[:, ip] - verts[:, i]
       le  = norm(e_v)
       n_v = e_v / le
       n_v = [n_v[2], -n_v[1]]
       for α = 1:2
            D[i,α] = integrate(ref1d,
            t -> dot(gmα[α+1]( verts[:, i] + t*e_v ), n_v))*le
            # D[i,α] = dot(gmα[α+1](0.), n_v)*le # if you wish the midpoint rule, comment out
       end
    end
    # stability constant
    σ    = 1.
    F    = Matrix(1.0I, ndof, ndof) - D*Π
    MASS = Π'*(Παβ*Π) + σ*(F' * F)

    # 3rd:  both parts have been computed
    return DIV, MASS
end
