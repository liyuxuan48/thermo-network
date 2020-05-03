module Systems

export ImpulseSystem,ShockTubeSystem,UUtoEverything

using ..Tools
"""
impulsesystem is a struct containing

    ρ = density
    G = mass flow rate
    P = perimeter of the cross section of the wall
    Ac = Area of the cross section of the wall
    qw = heat flux of the wall when t>0
    Δt = time interval
    Δz = space interval
"""

struct ImpulseSystem
    ρ
    G
    P
    Ac
    qw
    Δt
    Δz
end

"""
shocktubesystem is a struct containing

    gamma = gamma
    Δt = time interval
    Δz = space interval
"""

struct ShockTubeSystem
    gamma
    Δt
    Δz
end

"""
UUtoEverything is a struct containing all the thermodyamic values and flow velocity

    u = flow velocity
    ρ = density
    m = mass flow rate
    p = pressure
    e = total energy per mass
    ϵ = internal energy per mass
    h = enthalpy per mass

"""

struct UUtoEverything
    u
    ρ
    m
    p
    e
    ϵ
    h
    gamma
    c
end

function UUtoEverything(uu,gamma)

    u = uu[2,:]./uu[1,:]
    ρ = uu[1,:]
    m = u.*ρ
    p = uutourhop(uu,gamma)[3,:]
    e = uu[3,:]./ρ
    ϵ = p./ρ./(gamma-1)
    h = ϵ.+p./ρ
    c = sqrt.(gamma.*p./ρ)
return UUtoEverything(u,ρ,m,p,e,ϵ,h,gamma,c)

end

end
