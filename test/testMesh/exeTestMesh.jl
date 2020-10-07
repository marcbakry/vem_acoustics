include("../../mesh.jl")
include("../../utils.jl")
include("../../connectique.jl")
include("../../visu.jl")

function main()
    println("#### Validation maillage ####")
    # toy problem
    vtx = collect(transpose(
        Float64[1.0000    1.0000;
        0.9986    0.5823;
        0.8514    0.7718;
        0.8245    1.0000;
        0.7638    0.6373;
        0.7967    0.5113;
        0.5610    0.2393;
        1.0000    0.2083;
        1.0000   -0.0000;
        0.5663    0.9447;
        0.5567    1.0000;
        0.0000    1.0000;
        0.4946    0.2228;
        0.4691    0.0000;
        0.0000    0.5989;
       -0.0000    0.0000]
    ))
    elt = [
        Int64[15, 13, 10, 11, 12],
        Int64[3, 2, 1, 4],
        Int64[5, 3, 4, 11, 10],
        Int64[7, 8, 2, 6],
        Int64[13, 7, 6, 5, 10],
        Int64[16, 14, 13, 15],
        Int64[6, 2, 3, 5],
        Int64[14, 9, 8, 7, 13]
    ];
    # compute neighbours
    e2neigh, _ = getNeighbours(vtx, elt)
    nedg    = getNbEdges(e2neigh)
    # compute free nodes and edges 
    e2edg, freeedg, freenod = getMeshEdges(vtx, elt, nedg, e2neigh)
    #
    data = collect(range(0, stop=1., length=length(elt)));
    plotSolutionElt(mesh(vtx, elt, e2edg, freeedg, freenod, zeros(Int64, 2), nedg), data)
    nothing
end

main()