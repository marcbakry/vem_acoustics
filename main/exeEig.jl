#=
This code implements the VEM of order 0 for the computation of the resonance
modes of the Helmholtz equation. The method may be found in

"A virtual element method for the acoustic vibration problem"
Lourenço Beirão da Veiga, David Mora, Gonzalo Rivera, Rodolfo Rodríguez
https://arxiv.org/abs/1601.04316
=#

include("../utils.jl")  # some useful functions
include("../mesh.jl")   # defines a very general 2d mesh structure
include("../connectique.jl") # builds connectivity for the lower order method
include("../buildmatrices.jl") # builds the VEM matrices
include("../visu.jl") # visualization tools
include("../melem_rt0.jl") # elementary matrices with Raviart-Thomas of order 0

using LinearAlgebra  # package centralizing all linear algebra tools
using Arpack         # for the computation of the eigen values/vectors

# compute the pressure from u
function pressure(msh::mesh, u::Array{Float64, 1}, conn::Connectivity)
    p = zeros(Float64, length(msh.elt))
    for e in eachindex(msh.elt)
        E = polygonArea(msh.vtx[:,msh.elt[e]])
        for ee in conn.e2dof[e]
            sgnI = sign(ee)
            absI = abs(ee)
            p[e] += sgnI*u[absI]
        end
        p[e] /= -E
    end
    return p
end

function main()
    # read mesh
    println(" -> Reading mesh")
    # msh = readMesh("../meshes/voronoi.txt")
    # msh = readMesh("../meshes/squares.txt")
    # msh = readMesh("../meshes/triangles.txt")
    msh = readMesh("../meshes/non-convex.txt")
    # msh = readMesh("../meshes/smoothed-voronoi.txt")
    msh.vtx[2, :] .*= 1.1

    # connectivity
    println(" -> Building connectivity")
    conn = connectivity_k0(msh)

    # matrices
    println(" -> Assembling the matrices")
    A, B = buildmatrices_k0(msh, conn)
    # comment out below if you want to assemble with RT0; will only give a
    # correct solution for the mesh "triangles.txt"
    # A, B = buildRT0matrices(msh, conn)

    # apply boundary conditions
    mask = findall(msh.freeedg .== false)
    A = A[mask, mask]; # I vectorize because I am too lazy to write a loop
    B = B[mask, mask];

    # EV problem
    println(" -> Solving the generalized EV problem")
    Λ = eigen(A, B) # we get a structure with the eigenvalues and eigenvectors

    # plot eigenvectors
    it  = 0
    # index of the non-zero EV you wish to plot
    ind = 15
    for λ in sort(Λ.values)
        it += 1
        if λ > 1.05 # we only seek the non-zero eigenvalues
            uloc = Λ.vectors[:, it+ind]
            u    = zeros(Float64, conn.ndof)
            u[mask] .= uloc
            p    = pressure(msh, u, conn)
            plotSolutionElt(msh, p,
                myTitle="Pressure field - \$\\lambda=$(Λ.values[it+ind]-1)\$")
            #
            savefig("./EV_nonConvex.png")
            break
        end
    end
    println(Λ.values[it:(it+10)] .- 1)
end


main()
