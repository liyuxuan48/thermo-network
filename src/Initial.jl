module Initial

export initial_sod,initial_static,initial_uniform,initial_one_wave
"""
    this function's inputs are nx, xlim, and gamma

        nx is the number of grid points

        xlim is the limits of the domain

        gamma is the heat capacity ratio

    this function's outputs are ux and uu

        ux has 1 row
        the 1st row is x (x cooridnates)

        uu has 3 rows
        the 1st row is u (velocity)
        the 2nd row is m (mass flow rate)
        the 3rd row is ̂E (total energy per volume)

    this function initialize the Sod shock tube problem (Steger and Warming, 1981)
"""

function initial_sod(nx::Int64, xlim::Any, gamma::Float64)

    Lx = xlim[2]-xlim[1];
    Δx = Lx/(nx-1)

    # Initialization of ux and urhop
    ux = xlim[1]:Δx:xlim[2]
    urhop = Array{Float64,2}(UndefInitializer(), 3,length(ux))

    for i in 1:1:nx
    if (ux[i] < (xlim[2]-xlim[1])/2) urhop[:,i] = [0,1,1]
    else urhop[:,i] = [0,0.125,0.1] end
    end

    #From urhop make uu (u rho*u e)
    uu = Array{Float64,2}(UndefInitializer(), 3,length(ux))
    uu[1,:] = urhop[2,:];
    uu[2,:] = urhop[2,:].*urhop[1,:];
    uu[3,:] = urhop[3,:]/(gamma-1) + 0.5urhop[2,:].*urhop[1,:].^2;

    return ux,uu

end

"""
    this function's inputs are nx, xlim, and gamma

        nx is the number of grid points

        xlim is the limits of the domain

        gamma is the heat capacity ratio

    this function's outputs are ux and uu

        ux has 1 row
        the 1st row is x (x cooridnates)

        uu has 3 rows
        the 1st row is u (velocity)
        the 2nd row is m (mass flow rate)
        the 3rd row is ̂E (total energy per volume)

    this function initialize the Sod shock tube problem (Steger and Warming, 1981)
"""
function initial_static(nx::Int64, xlim::Any, gamma::Float64)

    Lx = xlim[2]-xlim[1];
    Δx = Lx/(nx-1)

    # Initialization of ux and urhop
    ux = xlim[1]:Δx:xlim[2]
    urhop = Array{Float64,2}(UndefInitializer(), 3,length(ux))

    for i in 1:1:nx
    if (ux[i] < (xlim[2]-xlim[1])/2) urhop[:,i] = [0,0.125,1]
    else urhop[:,i] = [0,0.125,1] end
    end

    #From urhop make uu (u rho*u e)
    uu = Array{Float64,2}(UndefInitializer(), 3,length(ux))
    uu[1,:] = urhop[2,:];
    uu[2,:] = urhop[2,:].*urhop[1,:];
    uu[3,:] = urhop[3,:]/(gamma-1) + 0.5urhop[2,:].*urhop[1,:].^2;

    return ux,uu

end

"""
    this function's inputs are nx, xlim, and gamma

        nx is the number of grid points

        xlim is the limits of the domain

        gamma is the heat capacity ratio

    this function's outputs are ux and uu

        ux has 1 row
        the 1st row is x (x cooridnates)

        uu has 3 rows
        the 1st row is u (velocity)
        the 2nd row is m (mass flow rate)
        the 3rd row is ̂E (total energy per volume)

    this function initialize the Sod shock tube problem (Steger and Warming, 1981)
"""
function initial_uniform(nx::Int64, xlim::Any, gamma::Float64)

    Lx = xlim[2]-xlim[1];
    Δx = Lx/(nx-1)

    # Initialization of ux and urhop
    ux = xlim[1]:Δx:xlim[2]
    urhop = Array{Float64,2}(UndefInitializer(), 3,length(ux))

    for i in 1:1:nx
    if (ux[i] < (xlim[2]-xlim[1])/2) urhop[:,i] = [1,0.125,1]
    else urhop[:,i] = [1,0.125,1] end
    end

    #From urhop make uu (u rho*u e)
    uu = Array{Float64,2}(UndefInitializer(), 3,length(ux))
    uu[1,:] = urhop[2,:];
    uu[2,:] = urhop[2,:].*urhop[1,:];
    uu[3,:] = urhop[3,:]/(gamma-1) + 0.5urhop[2,:].*urhop[1,:].^2;

    return ux,uu

end

end
