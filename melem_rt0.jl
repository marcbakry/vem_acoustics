#= to check : Raviart-Thomas assembling =#

function getBaseRT0(E::Float64, xE::Array{Float64, 2})
    return [
        x -> (x - xE[:, 3])/(2*E),
        x -> (x - xE[:, 1])/(2*E),
        x -> (x - xE[:, 2])/(2*E)
    ]
end

function melemRT0(vtx::Array{Float64, 2})
    @assert size(vtx, 2) == 3
    # initialize
    E = polygonArea(vtx)
    # div part
    DIV = ones(Float64, 3, 3) / E
    # mass part
    RT0 = getBaseRT0(E, vtx)
    MASS = zeros(Float64, 3, 3)
    quad = quadGL2d(3, vtx)
    for j=1:3
        for i=1:3
            MASS[i,j] = integrate(quad, x -> dot(RT0[i](x),RT0[j](x)))
        end
    end
    return DIV, MASS
end


function buildRT0matrices(msh::mesh, conn::Connectivity)
    e2dof = conn.e2dof
    ndof  = conn.ndof
    # 1st: preallocate
    A = zeros(Float64, ndof, ndof)
    B = zeros(Float64, ndof, ndof)
    # 2nd: main loop
    for e=1:length(msh.elt)
        @inbounds verts  = msh.vtx[:, msh.elt[e]]
        AE, BE = melemRT0(verts)
        nloc = length(e2dof[e])
        for i = 1:nloc
            @inbounds absI = abs(e2dof[e][i])
            @inbounds sgnI = sign(e2dof[e][i])
            for j=1:nloc
                @inbounds absJ = abs(e2dof[e][j])
                @inbounds sgnJ = sign(e2dof[e][j])
                @inbounds A[absI, absJ] += sgnI*sgnJ*(AE[i,j] + BE[i,j])
                @inbounds B[absI, absJ] += sgnI*sgnJ*BE[i,j]
                #
            end
        end
    end
    return A, B
end
