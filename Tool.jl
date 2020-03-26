"""
    this function's inputs are uu and gamma

        uu has 3 rows
        the 1st row is u (velocity)
        the 2nd row is m (mass flow rate)
        the 3rd row is e (total energy per mass)

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
        the 3rd row is e (total energy per mass)

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
