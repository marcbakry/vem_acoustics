#---------------------------------------------#
# Geometrical transformations for the segment #
#---------------------------------------------#

function getSegTransform(A::Array{Float64, 1}, B::Array{Float64,1})
    #= Returns (u0, u) the transform data from [-1, 1] to a 2d segment [A, B] =#
    return (A + B)/2, (B-A)/2
end

 @inline function segRefToSeg2d(xref::Float64, A::Array{Float64,1},
    B::Array{Float64, 1})
    #= Transformation from the 1d [-1, 1] reference segment to a 2-dimensional
        segment [A, B] where A and B are the nodes =#
    return A + (xref + 1)/2*(B-A)
end

@inline function segRefToSeg2d(xref::Float64, u0::Array{FLoat64, 1},
    u::Array{Float64})
    #= Affine transformation ans = u0 + xref*u where u0 and u are obtained from
        getSegTransform() =#
    return u0 + xref*u
end

#----------------------------------------------#
# Geometrical transformations for the triangle #
#----------------------------------------------#
function getTriTransform(T::Array{Float64,2})
    #= Returns the matrices required for the computation of the transformation
        from the triangle {(0,0),(1,0),(1,0)} to {A, B, C} where A, B and C are
        stored columnwise in T =#
    return T[:, 1], [T[1,2]-T[1,1] T[1,3]-T[1,1];T[2,2]-T[2,1] T[2,3]-T[2,1]]
end

@inline function triRefToAny(xref::Array{Float64,1}, T::Array{Float64, 2})
    #= Transformation from the reference triangle {(0,0),(1,0),(1,0)} to
        T={A, B, C} =#
    return T[:, 1] + [T[1,2]-T[1,1] T[1,3]-T[1,1];T[2,2]-T[2,1] T[2,3]-T[2,1]]*xref
end

@inline function triRefToAny(xref::Array{Float64}, u0::Array{Float64,1},
    u::Array{Float64})
    #= Transformation from the reference triangle {(0,0),(1,0),(1,0)} to
        T={A, B, C} using the matrices returned by getTriTransform() =#
    return u0 + u*xref
end
