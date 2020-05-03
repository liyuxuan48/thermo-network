module BoundaryCondition

export set_outlet_nonreflect_boundary!

using..Systems


    """
        get the characteristic wave amplitudes L
    """

    function get_L_from_nonreflect(uu::Array,everythinginitial,Δx::Float64)

            # import variables
            gamma = everythinginitial.gamma


            uueverything = UUtoEverything(uu,gamma)

            u = uueverything.u
            ρ = uueverything.ρ
            c = uueverything.c
            p = uueverything.p

            # get λ1,λ2,λ3,λ4,λ5
            λ = Array{Float64,1}(UndefInitializer(), 5)
            λ[1] = u[end]-c[end]
            λ[2] = u[end]
            λ[3] = u[end]
            λ[4] = u[end]
            λ[5] = u[end]+c[end]

            # get L1,L2,L3,L4,L5
            L = Array{Float64,1}(UndefInitializer(), 5)
            L[1]=λ[1].*0
            L[2]=λ[2].*(c[end].^2 .* (ρ[end]-ρ[end-1])./Δx-(p[end]-p[end-1])./Δx)
            L[3]=λ[3].*0
            L[4]=λ[4].*0
            L[5]=λ[5].*((p[end]-p[end-1])./Δx+ρ[end].*c[end].*(u[end]-u[end-1])./Δx)
#             println("L=",L)

        return L
        end

        """
            get the d from L
        """

        function get_d_from_L(uu::Array,everythinginitial,L::Array)

        # import variables
        gamma = everythinginitial.gamma

        uueverything = UUtoEverything(uu,gamma)

        u = uueverything.u
        ρ = uueverything.ρ
        c = uueverything.c
        p = uueverything.p

        # get d1,d2,d3,d4,d5
        d = Array{Float64}(UndefInitializer(), 5)
        d[1] = 1 ./ (c[end].^2).*(L[2]+0.5.*(L[5]+L[1]))
        d[2] = 0.5 .* (L[5]+L[1])
        d[3] = 0.5 ./ρ[end]./c[end] .* (L[5]-L[1])
        d[4] = 0
        d[5] = 0

        return d
        end


    """
        still working on this
    """
    function set_outlet_nonreflect_boundary!(uu::Array,everythinginitial,Δx::Float64,Δt::Float64)

    L = get_L_from_nonreflect(uu,everythinginitial,Δx)
    d = get_d_from_L(uu,everythinginitial,L)

    gamma = everythinginitial.gamma
    h = everythinginitial.h

    uueverything = UUtoEverything(uu,gamma)

    u = uueverything.u
    ρ = uueverything.ρ

    uuend=Array{Float64,1}(UndefInitializer(), 3)
    uuend[1]=uu[1,end] + (-d[1]).*Δt
    uuend[2]=uu[2,end] + (-u[end].*d[1]-ρ[end].*d[3]+0).*Δt
    uuend[3]=uu[3,end] + (-0.5 .* u[end].*u[end].*d[1]-d[2]./(gamma-1) - ρ[end].*u[end].*d[3] + 0).*Δt

    return uuend
    end


# """
#     still working on this
# """
# """
#
#     function set_h_boundary!(uu::Array,everythinginitial)
#
#         gamma = everythinginitial.gamma
#         h = everythinginitial.h
#
#
#         uueverything = UUtoEverything(uu,gamma)
#
#         u = uueverything.u
#         ρ = uueverything.ρ
#
#         ϵ = h./gamma
#         Ehat = ρ.*ϵ + 0.5.*ρ.*u.*u
#
#         uunew = Array{Float64,2}(UndefInitializer(), 3,size(uu)[2])
#
#         uunew[1,:]=uu[1,:]
#         uunew[2,:]=uu[2,:]
#         uunew[3,:]=Ehat
#
#     return uunew[:,1]
#     end
# """
#
#
#     """
#     function setuuboundary!(uu::Array,everythinginitial::UUtoEverything)
#
#         gamma = everythinginitial.gamma
#         h = everythinginitial.h
#
#
#         uueverything = UUtoEverything(uu,gamma)
#
#         u = uueverything.u
#         ρ = uueverything.ρ
#
#         ϵ = h./gamma
#         Ehat = ρ.*ϵ + 0.5.*ρ.*u.*u
#
#         uunew = Array{Float64,2}(UndefInitializer(), 3,size(uu)[2])
#
#         uunew[1,:]=everythinginitial.ρ
#         uunew[2,:]=everythinginitial.m
#         uunew[3,:]=Ehat
#
#     return uunew[:,1]
#     end
#     """
end
