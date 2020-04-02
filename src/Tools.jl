module Tools

export uutourhoc,uutourhop,heaviside,set_h_boundary!

"""
    this function's inputs are uu and gamma

        uu has 3 rows
        the 1st row is u (velocity)
        the 2nd row is m (mass flow rate)
        the 3rd row is ̂E (total energy per volume)

        gamma is the heat capacity ratio

    this function's outputs is tourhoc

        tourhoc has 3 rows
        the 1st row is u (velocity)
        the 2nd row is ρ (density)
        the 3rd row is c (speed of sound)

    this function uses perfect gas state function to get tourhoc
"""
function uutourhoc(uu::Any, gamma::Float64)

    # tourhoc = Array{Float64,2}(UndefInitializer(),length(uu[:,1]),length(uu[1,:]))

    tourhoc = deepcopy( zeros(size(uu)))

    tourhoc[1,:] = uu[2,:]./uu[1,:]
    tourhoc[2,:] = uu[1,:]

    pp = ((uu[3,:]-0.5tourhoc[2,:].*tourhoc[1,:].^2)*(gamma-1))

    tourhoc[3,:] = sqrt.(gamma.*abs.(pp)./(tourhoc[2,:]))

    return tourhoc

end

"""
    this function's inputs are uu and gamma

        uu has 3 rows
        the 1st row is u (velocity)
        the 2nd row is m (mass flow rate)
        the 3rd row is ̂E (total energy per volume)

        gamma is the heat capacity ratio

    this function's outputs is tourhop

        tourhoc has 3 rows
        the 1st row is u (velocity)
        the 2nd row is ρ (density)
        the 3rd row is p (pressure)

    this function uses perfect gas state function to get tourhop
"""

function uutourhop(uu::Any, gamma::Float64)

    tourhop = deepcopy( zeros(size(uu)))

    tourhop[1,:] = uu[2,:]./uu[1,:]
    tourhop[2,:] = uu[1,:]

    tourhop[3,:] = ((uu[3,:]-0.5tourhop[2,:].*tourhop[1,:].^2)*(gamma-1))

    return tourhop

end

"""
this function's inputs is x

    x is the coordinate

this function's outputs is a Float value

    this function is the heaviside step function (half-maximum convention)
"""
heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))

end
