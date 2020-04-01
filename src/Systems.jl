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
