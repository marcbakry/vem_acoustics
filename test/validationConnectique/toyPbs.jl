#= Generate toy problems =#

function generateToyProblem1()
    vtx = collect(transpose(
        Float64[0 0; 2 0; 4 0; 4 2; 4 4; 2 4; 0 4; 0 2; 2 1; 3 2; 2 3; 1 2]))
    bnd = Int64[1,2,3,4,5,6,7,8]
    elt = [Int64[1,2,9,12,8],
        Int64[2,3,4,10,9],
        Int64[10,4,5,6,11],
        Int64[8,12,11,6,7],
        Int64[9,10,11,12]
    ]
    return vtx, elt, bnd
end


function generateToyProblem2()
    vtx = collect(transpose(Float64[0 0; 0.5 0; 0 1; 1 0; 2 0; 2 1; 1.5 1; 1.5 0.5; 1 1.5]))
    bnd = Int64[1,2,3,4,5,6,7,9]
    elt = [
        Int64[1,2,3],
        Int64[2,4,8,7,9,3],
        Int64[4,5,6,7,8]
    ]
    return vtx, elt, bnd
end

function correctDataToyProblem2()
    nod2e = [
        Int64[1], Int64[1, 2], Int64[1, 2], Int64[2,3], Int64[3], Int64[3],
        Int64[2, 3], Int64[2,3], Int64[2]
    ]
    e2neighs = [
        Int64[0, 2, 0], Int64[0, 3, 3, 0, 0, 1], Int64[0, 0, 0, 2, 2]
    ]
    nedg = Int64(11)
    e2edg = [
        Int64[1, 2, 3], Int64[4, 5, 6, 7, 8, 2], Int64[9, 10, 11, 6, 5]
    ]
    freeedg = [true, false, true, true, false, false, true, true, true, true, true]
    freenod = [true, true, true, true, true, true, true, false, true]
    return nod2e, e2neighs, nedg, e2edg, freeedg, freenod
end

function connectivityToyProblem2()
    conn = [
        Int64[1,2,3], Int64[4,5,6,7,8,-2], Int64[9,10,11,-6,-5]
    ]
    return conn
end
