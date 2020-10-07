#
include("mesh.jl")

mutable struct Connectivity
    ndof::Int64
    e2dof::Array{Array{Int64,1},1}
    bnddofs::Array{Int64, 1}
end

function connectivity_k0(msh::mesh)
    #= Compute the connectivity table for the cas k=0 =#
    e2edg = msh.edg
    e2nod = msh.elt
    nedg  = msh.nedg
    nelt = length(e2nod)
    conn = [zeros(Int64, length(e2edg[i])) for i=1:nelt]
    edgflagged = zeros(Bool, nedg)
    for e=1:nelt
        edgloc  = e2edg[e]
        nedgloc = length(edgloc)
        for ee=1:nedgloc
            if edgflagged[edgloc[ee]] != true
                conn[e][ee] = edgloc[ee]
                edgflagged[edgloc[ee]] = true
            else
                conn[e][ee] = -edgloc[ee]
            end
        end
    end
    # dofs on the boundary
    nbnddof = sum(msh.freeedg)
    bnddofs = zeros(Int64, nbnddof)
    it = 0
    for e in eachindex(msh.freeedg)
        if msh.freeedg[e]
            it += 1
            bnddofs[it] = e
        end
    end
    #
    return return Connectivity(nedg, conn, bnddofs)
end
