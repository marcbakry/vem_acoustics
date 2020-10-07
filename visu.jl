#=
My very basic visualization toolkit for the VEM using the Plots module.
=#

using Plots
pyplot() # select PyPlot background

#= The next three functions are used to produce the Jet colormap. =#
function interpolate(v::Float64, y0::Float64, x0::Float64, y1::Float64, x1::Float64)
    return (v - x0)*(y1 - y0)/(x1 - x0) + y0
end

function base(v::Float64)
    if v <= -0.75
        return 0.
    elseif v <= -0.25
        return interpolate(v, 0., -0.75, 1., -0.25)
    elseif v <= 0.25
        return 1.
    elseif v <= 0.75
        return interpolate(v, 1., 0.25, 0., 0.75)
    else
        return 0.
    end
end

red(gs::Float64)   = base(gs - 0.5)
green(gs::Float64) = base(gs)
blue(gs::Float64)  = base(gs + 0.5)

function plotSolutionNodal(msh::mesh, u::Array{Float64, 1}; myTitle::String="",
    toFile::Bool=false)
    #= Plots nodal values =#
    println("Printing solution")
    # set polygons
    nvtx = size(msh.vtx, 2)
    nelt = length(msh.elt)
    moyU = zeros(Float64, nelt)
    # set polygons
    el = msh.elt[1]
    xl = collect(transpose(msh.vtx[:,el]))
    println
    pols = [xl; NaN NaN]
    ul   = u[el]
    moyU[1] = sum(ul)/length(ul)
    for i=2:nelt
        @inbounds el      = msh.elt[i]
        @inbounds xl      = collect(transpose(msh.vtx[:, el]))
        tmp = [xl; NaN NaN]
        pols = [pols; tmp]
        @inbounds ul      = u[el]
        @inbounds moyU[i] = sum(ul)/length(ul)
    end
    # set colors
    maxU = maximum(moyU)
    minU = minimum(moyU)
    delU = maxU - minU
    gs   = -1 + 2/delU*(moyU[1] - minU)
    colU = [RGBA(red(gs), green(gs), blue(gs), 1)]
    for i=2:nelt
        @inbounds gs = -1 + 2/delU*(moyU[i] - minU)
        push!(colU, RGBA(red(gs), green(gs), blue(gs), 1))
    end
    # plot
    fig = Plots.plot(pols[:,1], pols[:,2], seriestype=:shape, c=colU, show=true,
    label=nothing, linewidth=0.2, linecolor=RGBA(255/255,255/255,255/255,1),
    aspect_ratio=:equal, title=myTitle)
    if toFile
        png(fig, "myVEM")
    end
    display(fig)
    nothing
end

function plotSolutionElt(msh::mesh, ul::Array{Float64, 1}; myTitle::String="",
    toFile::Bool=false)
    #= Plots a function constant per element, for example the pressure field. =#
    println("Printing solution")
    # set polygons
    nvtx = size(msh.vtx, 2)
    nelt = length(msh.elt)
    moyU = zeros(Float64, nelt)
    # set polygons
    el = msh.elt[1]
    xl = collect(transpose(msh.vtx[:,el]))
    println
    pols = [xl; NaN NaN]
    for i=2:nelt
        @inbounds el      = msh.elt[i]
        @inbounds xl      = collect(transpose(msh.vtx[:, el]))
        tmp = [xl; NaN NaN]
        pols = [pols; tmp]
    end
    # set colors
    maxU = maximum(ul)
    minU = minimum(ul)
    delU = maxU - minU
    gs   = -1 + 2/delU*(ul[1] - minU)
    colU = [RGBA(red(gs), green(gs), blue(gs), 1)]
    for i=2:nelt
        @inbounds gs = -1 + 2/delU*(ul[i] - minU)
        push!(colU, RGBA(red(gs), green(gs), blue(gs), 1))
    end
    #
    fig = Plots.plot(pols[:,1], pols[:,2], seriestype=:shape, c=colU, show=true,
    label=nothing, linewidth=0.2, linecolor=RGBA(255/255,255/255,255/255,1),
    aspect_ratio=:equal, title=myTitle)
    if toFile
        png(fig, "myVEM")
    end
    display(fig)
    nothing
end
