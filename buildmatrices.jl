#
include("mesh.jl")
include("melem.jl")

function buildmatrices_k0(msh::mesh, conn::Connectivity)
    #=
    Builds the matrix in full format (because I wanted to use the exact EV
    solver provided by Arpack).
    =#
    e2dof = conn.e2dof
    ndof  = conn.ndof
    # 1st: preallocate
    A = zeros(Float64, ndof, ndof)
    B = zeros(Float64, ndof, ndof)
    # 2nd: main loop
    for e=1:length(msh.elt)
        # @inbounds removes the check on the bounds of the array
        @inbounds verts  = msh.vtx[:, msh.elt[e]] # get the vertices
        # Choose if you want to compute the elementary matrix by projecting on
        # ∇P1 or on ∇P2 (for example for comparison with the RT0 :) )
        AE, BE = melemk0P1(verts)
        # AE, BE = melemk0P2(verts)
        nloc = length(e2dof[e])
        for i = 1:nloc
            @inbounds tmpI = e2dof[e][i]
            absI, sgnI = abs(tmpI), sign(tmpI)
            for j=1:nloc
                @inbounds tmpJ = e2dof[e][j]
                absJ, sgnJ = abs(tmpJ), sign(tmpJ)
                @inbounds A[absI, absJ] += sgnI*sgnJ*(AE[i,j] + BE[i,j])
                @inbounds B[absI, absJ] += sgnI*sgnJ*BE[i,j]
            end
        end
    end
    return A, B
end
