#=
THis file contains all my mesh-manadgement methods.
=#

# mesh related functions
include("utils.jl")

# MESH MANAGEMENT
mutable struct mesh
    vtx::Array{Float64, 2} # vertices
    elt::Array{Array{Int64, 1}, 1} # elt → vtx
    edg::Array{Array{Int64, 1}, 1} # elt → edg
    freeedg::Array{Bool, 1} # elt → free edge
    freenod::Array{Bool, 1} # elt → free node
    bnd::Array{Int64, 1} # list of boundary nodes
    nedg::Int64 # number of edges
end

function searchCommonElt(vtx2elt::Array{Array{Int64,1},1}, n1::Int64,
    n2::Int64, egs::Int64)
    #= Search for the common element. This function may be shortened by making
        better use of the language… =#
    nv2e1 = length(vtx2elt[n1])
    nv2e2 = length(vtx2elt[n2])
    egs0 = 0
    egs1 = 0
    egs2 = 0
    ind0 = 0
    j    = 1
    while ind0 == 0 && j <= nv2e1
        egs1 = vtx2elt[n1][j] # j-th elt of node n1
        if egs1 != egs
            ind = 0; k = 1
            while ind == 0 && k <= nv2e2
                egs2 = vtx2elt[n2][k]  # the candidate
                if egs1 == egs2
                    ind = 1
                else
                    k += 1
                end
            end
            if ind == 1
                ind0 = 1
            else
                j += 1
            end
        else
            j += 1
        end
    end
    #
    if ind0 == 1
        if egs1 == egs2 && egs1 != egs
            egs0 = egs1
        else
            egs0 = 0
        end
    else
        egs0 = 0
    end
    #
    return egs0
end

function getNeighbours(vtx::Array{Float64, 2}, elt::Array{Array{Int64,1},1})
    #= This function returns the elements on each side of an edge. The value 0
        means that there is only one element ≡ the edge is free =#
    nvtx = size(vtx, 2)
    nelt = length(elt)
    # 1st: for each node, determine all neighbour elements
    vtx2elt = [Int64[] for i=1:nvtx]
    for ie=1:nelt # loop over all elements
        for ix in elt[ie] # loop over all nodes of each element
            push!(vtx2elt[ix], ie) # add the element to each of these nodes
        end
    end
    # 2nd: set the correspondance between the edges and the elements
    e2neighs = [Int64[] for i=1:nelt]
    for ie=1:nelt
        neltie = length(elt[ie])
        for k=1:neltie
            kp  = mod_wrap(k+1, neltie)
            nk  = elt[ie][k]
            nkp = elt[ie][kp]
            egs0 = searchCommonElt(vtx2elt, nk, nkp, ie)
            push!(e2neighs[ie], egs0)
        end
    end
    #
    return e2neighs, vtx2elt
end

function getNbEdges(e2neighs::Array{Array{Int64,1},1})
    #= Compute the total number of edges =#
    nelt = length(e2neighs)
    flag = zeros(Bool, nelt)
    nedg = Int64(0)
    for e=1:nelt
        neloc = length(e2neighs[e])
        for eloc=1:neloc
            tmp = e2neighs[e][eloc]
            if tmp != 0
                if flag[tmp] != true
                    nedg += 1
                end
            else
                nedg += 1
            end
        end
        flag[e] = true
    end
    #
    return nedg
end

function getMeshEdges(vtx::Array{Float64, 2}, elt::Array{Array{Int64, 1}, 1},
    nedg::Int64, e2neighs::Array{Array{Int64,1},1})
    #= This function assigns numbers to the edges of the mesh and return the
        correspondance between elements and edges with a sign =#
    nnod = size(vtx, 2)
    nelt = length(elt)
    #
    freeedg = zeros(Bool, nedg)
    freenod = zeros(Bool, nnod)
    nedg = 0 # reinitialize
    e2edg = [Int64[] for i=1:nelt]
    for e=1:nelt
        nodloc   = elt[e]
        e2edg[e] = zeros(Int64,length(nodloc))
    end
    #
    for e=1:nelt
        nodloc   = elt[e]
        xt       = vtx[:, nodloc]
        for eloc = 1:length(nodloc)
            if e2edg[e][eloc] != 0
                continue
            end
            if e2neighs[e][eloc] != 0
                elocp = mod_wrap(eloc+1, length(nodloc))
                xmt   = (xt[:,eloc] + xt[:, elocp])/2
                nodlocneigh = elt[e2neighs[e][eloc]]
                xl    = vtx[:, nodlocneigh]
                for u=1:length(nodlocneigh)
                    up  = mod_wrap(u+1, length(nodlocneigh))
                    xml = (xl[:,u] + xl[:, up])/2
                    if norm(xmt - xml) < 1e-8
                        nedg += 1
                        e2edg[e][eloc] = nedg
                        e2edg[e2neighs[e][eloc]][u] = nedg
                        break
                    end
                end
            else
                nedg += 1
                e2edg[e][eloc] = nedg
                freeedg[nedg]  = true
            end
        end
    end
    # free nodes
    for e=1:nelt
        neloc = length(elt[e])
        for eloc=1:neloc
            if freeedg[e2edg[e][eloc]]
                freenod[elt[e][eloc]] = true
                freenod[elt[e][mod_wrap(eloc+1, neloc)]] = true
            end
        end
    end
    return e2edg, freeedg, freenod
end

function readMesh(fname::String)
    println("Reading mesh file: $fname")
    vtx = zeros(Float64, 1, 1)
    elt = Array{Int64, 1}[];
    bnd = zeros(Int64, 1)
    open(fname, "r") do fid
        ln = readline(fid)
        println(" -> Reading vertices ...")
        # vertices
        nvtx = parse(Int64, ln)
        vtx = zeros(Float64, 2, nvtx)
        for i=1:nvtx
            ln = readline(fid);
            vtx[:, i] .= parse.(Float64, split(ln, " ", keepempty=false))
        end
        #elements
        println(" -> Reading elements ...")
        ln = readline(fid)
        nelt = parse(Int64, ln)
        for i=1:nelt
            ln = readline(fid)
            ln_s = split(ln, " ", keepempty=false)
            push!(elt, parse.(Int64, ln_s))
        end
        # boundary
        println(" -> Reading boundary nodes ...")
        ln = readline(fid)
        nbnd = parse(Int64, ln)
        bnd = zeros(Int64, nbnd)
        for i=1:nbnd
            ln = readline(fid)
            bnd[i] = parse(Int64, ln)
        end
    end
    #
    e2neigh, _ = getNeighbours(vtx, elt)
    nedg    = getNbEdges(e2neigh)
    e2edg, freeedg, freenod = getMeshEdges(vtx, elt, nedg, e2neigh)
    return mesh(vtx, elt, e2edg, freeedg, freenod, bnd, nedg)
end

# compute connectivity table
function connectivityTableOrde0(msh::mesh)
    #= This function returns the connectivity table corresponding to a problem
        of order k = 0. =#
end
