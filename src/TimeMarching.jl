module TimeMarching

export  stegerwarmingrk1!,upwindrk1!

using ..Tools
using ..FiniteDifference

"""
    this function's inputs are t and uu

        t is the current time

        uu has 3 rows
        the 1st row is u (velocity)
        the 2nd row is m (mass flow rate)
        the 3rd row is e (total energy per mass)

    this function's outputs are new t and uunew

        t is the current time

        uunew has the same structure as uu

    this function uses flux splitting scheme (Steger and Warming, 1981)
    and 1st order Euler time marching

    The current function is only good for u<=c (speed of sound)
"""

function stegerwarmingrk1!(t,uu,righthand,shocktubesystem)

    gamma=shocktubesystem.gamma
    Δt=shocktubesystem.Δt
    Δx=shocktubesystem.Δz


    #u<=c
    urhoc = uutourhoc(uu,gamma)

    u = urhoc[1,:]
    rho = urhoc[2,:]
    c = urhoc[3,:]

    #initialization of uunew
    uunew =  deepcopy( zeros(size(uu)))

    #initialization of f
    fplus = deepcopy( zeros(size(uu)))
    fminus = deepcopy( zeros(size(uu)))


    #get fplus and fminus
    fplus[1,:] = (rho./(2gamma)).*(2gamma.*u+c-u)
    fplus[2,:] = (rho./(2gamma)).*(2(gamma-1).*u.^2+(u+c).^2)
    fplus[3,:] = (rho./(2gamma)).*((gamma-1).*u.^2+(((u+c).^3)/2)+((3-gamma).*(u+c).*c.^2)./(2(gamma-1)))

    fminus[1,:] = (rho./(2gamma)).*(u-c)
    fminus[2,:] = (rho./(2gamma)).*(u-c).^2
    fminus[3,:] = (rho./(2gamma)).*((((u-c).^3)/2)+((3-gamma).*(u-c).*c.^2)./(2(gamma-1)))

    for i in 2:1:length(uu[1,:])-1 # shrink a grid
        uunew[:,i] = (-(fplus[:,i]-fplus[:,i-1])/Δx)*Δt + (-(fminus[:,i+1]-fminus[:,i])/Δx)*Δt .+ righthand*Δt + uu[:,i]
    end

    # deal with the side
    uunew[:,length(uunew[1,:])] = uunew[:,length(uunew[1,:])-1]
    uunew[:,1] = uunew[:,2]

    t = t + Δt
    return t,uunew
end

"""
    this function's inputs are t and w and impulsesystem

        t is the current time

        w has 1 row
        the 1st row is h (enthalpy)

        impulsesystem is an ImpulseSystem struct

    this function's outputs are new t and out

        t is the current time

        out has the same structure as w

    this function uses upwind scheme and 1st order Euler time marching

    The current function is only good for information tavelling to the right
"""

function upwindrk1!(t::Any,w::Any,impulsesystem::Any)

    ρ = impulsesystem.ρ
    G = impulsesystem.G
    P = impulsesystem.P
    Ac = impulsesystem.Ac
    qw = impulsesystem.qw
    Δt = impulsesystem.Δt
    Δz = impulsesystem.Δz

    out = deepcopy(w) .+ (-G/ρ.*upwind(w,Δz) .+ P*qw/Ac/ρ).*Δt

    t = deepcopy(t) .+ Δt
    return t,out
end

end
