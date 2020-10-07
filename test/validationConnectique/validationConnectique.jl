include("../../mesh.jl")
include("../../utils.jl")
include("../../connectique.jl")
# include("../../visu.jl")

include("toyPbs.jl")

function validateToyProblem(nod2e::Array{Array{Int64, 1}, 1},
    e2neighs::Array{Array{Int64, 1}, 1}, nedg::Int64, e2edg::Array{Array{Int64, 1}, 1},
    freeedg::Array{Bool, 1}, freenod::Array{Bool, 1}, nod2e_ref::Array{Array{Int64, 1}, 1},
    e2neighs_ref::Array{Array{Int64, 1}, 1}, nedg_ref::Int64,
    e2edg_ref::Array{Array{Int64, 1}, 1}, freeedg_ref::Array{Bool, 1}, freenod_ref::Array{Bool, 1})
    # checking nod2e
    println("---------------------------------------------------")
    @assert length(nod2e)  == length(nod2e_ref)
    for n in eachindex(nod2e)
        println("Checking node: $n")
        @assert length(nod2e[n]) == length(nod2e_ref[n])
        for ne in eachindex(nod2e[n])
            println("  -> $(nod2e[n][ne])    -> $(nod2e_ref[n][ne])")
        end
    end
    # checking e2neighs
    println("---------------------------------------------------")
    @assert length(e2neighs) == length(e2neighs_ref)
    for e in eachindex(e2neighs)
        println("Checking element: $e")
        @assert length(e2neighs[e]) == length(e2neighs_ref[e])
        for ee in eachindex(e2neighs[e])
            println("  -> $(e2neighs[e][ee])    -> $(e2neighs_ref[e][ee])")
        end
    end
    # checking nedg
    println("---------------------------------------------------")
    println("Checking nedg: $nedg     $nedg_ref")
    # checking
    println("---------------------------------------------------")
    @assert length(e2edg) == length(e2edg_ref)
    for e in eachindex(e2edg)
        println("Checking element: $e")
        @assert length(e2edg[e]) == length(e2edg_ref[e])
        for ee in eachindex(e2edg[e])
            println("  -> $(e2edg[e][ee])    -> $(e2edg_ref[e][ee])")
        end
    end
    # check free edges
    println("---------------------------------------------------")
    println("Checking free edges")
    for e in eachindex(freeedg)
        println("  -> $(freeedg[e])    -> $(freeedg_ref[e])")
    end
    println("---------------------------------------------------")
    println("Checking free nodes")
    for n in eachindex(freenod)
        println("  -> $(freenod[n])    -> $(freenod_ref[n])")
    end
    nothing
end

function mainValidation2()
    vtx, elt, bnd = generateToyProblem2()
    # get mesh data
    e2neigh, nod2e = getNeighbours(vtx, elt)
    nedg           = getNbEdges(e2neigh)
    e2edg, freeedg, freenod = getMeshEdges(vtx, elt, nedg, e2neigh)
    # get correct mesh data
    nod2e_ref, e2neigh_ref, nedg_ref, e2edg_ref, freeedg_ref, freenod_ref =
        correctDataToyProblem2()
    # now check
    validateToyProblem(nod2e, e2neigh, nedg, e2edg, freeedg, freenod,
        nod2e_ref, e2neigh_ref, nedg_ref, e2edg_ref, freeedg_ref, freenod_ref)
    # connectivity
    msh = mesh(vtx, elt, e2edg, freeedg, freenod, bnd, nedg)
    conn = connectivity_k0(msh)
    conn_ref = connectivityToyProblem2()
    println("---------------------------------------------------")
    println("Checking connectivity")
    @assert length(conn.e2dof) == length(conn_ref)
    for e in eachindex(conn.e2dof)
        @assert length(conn.e2dof[e]) == length(conn_ref[e])
        for nd in eachindex(conn.e2dof[e])
            println("  -> $(conn.e2dof[e][nd])    -> $(conn_ref[e][nd])")
        end
    end
    #
    nothing
end

mainValidation2()
