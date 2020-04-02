module BoundaryCondition

export set_h_boundary!

using..Systems

    function set_h_boundary!(uu::Array,everythinginitial)

        gamma = everythinginitial.gamma
        h = everythinginitial.h


        uueverything = UUtoEverything(uu,gamma)

        u = uueverything.u
        ρ = uueverything.ρ

        ϵ = h./gamma
        e = ρ.*ϵ + 0.5.*ρ.*u.*u

        uunew = Array{Float64,2}(UndefInitializer(), 3,size(uu)[2])

        uunew[1,:]=uu[1,:]
        uunew[2,:]=uu[2,:]
        uunew[3,:]=e

    return uunew[:,1]
    end

    """
    function setuuboundary!(uu::Array,everythinginitial::UUtoEverything)

        gamma = everythinginitial.gamma
        h = everythinginitial.h


        uueverything = UUtoEverything(uu,gamma)

        u = uueverything.u
        ρ = uueverything.ρ

        ϵ = h./gamma
        e = ρ.*ϵ + 0.5.*ρ.*u.*u

        uunew = Array{Float64,2}(UndefInitializer(), 3,size(uu)[2])

        uunew[1,:]=everythinginitial.ρ
        uunew[2,:]=everythinginitial.m
        uunew[3,:]=e

    return uunew[:,1]
    end
    """
end
