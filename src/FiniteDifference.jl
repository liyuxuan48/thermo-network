module FiniteDifference

export upwind,downwind,central

"""
    this function's inputs are w and Δz

        w has 1 row
        the 1st row is h (enthalpy)

        Δz is the space interval

    this function's output is out

        out has the same structure as w

    this function uses upwind finite difference scheme to get the 1st order spacial derivative

    The current function is only good for information tavelling to the right
"""

function upwind(w::Array,Δz::Float64)
    out=zero(deepcopy(w))
    NX=length(w)
    @inbounds for x in 2:NX-1
        out[x] =  (w[x] - w[x-1])/Δz
    end

    out[1]=out[2]
    out[NX]=out[NX-1]

    return out
end

"""
    this function's inputs are w and Δz

        w has 1 row
        the 1st row is h (enthalpy)

        Δz is the space interval

    this function's output is out

        out has the same structure as w

    this function uses downwind finite difference scheme to get the 1st order spacial derivative

    The current function is only good for information tavelling to the left
"""
function downwind(w::Array,Δz::Float64)
    out=zero(deepcopy(w))
    NX=length(w)
    @inbounds for x in 2:NX-1
        out[x] =  (w[x+1] - w[x])/Δz
    end

    out[1]=out[2]
    out[NX]=out[NX-1]

    return out
end

"""
this function's inputs are w and Δz

    w has 1 row
    the 1st row is h (enthalpy)

    Δz is the space interval

this function's output is out

    out has the same structure as w

this function uses central finite difference scheme to get the 2st order spacial derivative
"""

function central(w::Array,Δz::Float64)
    out=zero(deepcopy(w))
    NX=length(w)
    @inbounds for x in 2:NX-1
        out[x] =  (w[x+1] - w[x-1])/(2Δz)
    end

    out[1]=out[2]
    out[NX]=out[NX-1]

    return out
end

end
