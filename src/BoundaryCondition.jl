module BoundaryCondition

export set_outlet_nonreflect_boundary!,set_outlet_costant_p_boundary!,
set_inlet_constant_h!,set_outlet_costant_h_boundary!,set_inlet_interface!,set_outlet_interface!,
get_L_from_interface,get_L_from_nonreflect,get_d_from_L_inflow,get_d_from_L_outflow

using..Systems

function get_L_from_interface(uu1::Array,uu2::Array,everythinginitial1,everythinginitial2,Δx::Float64)

            # import variables
            gamma = everythinginitial1.gamma


            uueverything1 = UUtoEverything(uu1,gamma)

            u1 = uueverything1.u
            ρ1 = uueverything1.ρ
            c1 = uueverything1.c
            p1 = uueverything1.p

            # get λ1,λ2,λ3,λ4,λ5
            λ1 = Array{Float64,1}(UndefInitializer(), 5)
            λ1[1] = u1[end]-c1[end]
            λ1[2] = u1[end]
            λ1[3] = u1[end]
            λ1[4] = u1[end]
            λ1[5] = u1[end]+c1[end]


            # get L2,L3,L4,L5 from upstream
            L = Array{Float64,1}(UndefInitializer(), 5)
            L[2]=λ1[2].*(c1[end].^2 .* (ρ1[end]-ρ1[end-1])./Δx-(p1[end]-p1[end-1])./Δx)
            L[3]=λ1[3].*0
            L[4]=λ1[4].*0
            L[5]=λ1[5].*((p1[end]-p1[end-1])./Δx+ρ1[end].*c1[end].*(u1[end]-u1[end-1])./Δx)


            uueverything2 = UUtoEverything(uu2,gamma)

            u2 = uueverything2.u
            ρ2 = uueverything2.ρ
            c2 = uueverything2.c
            p2 = uueverything2.p

            # get λ1,λ2,λ3,λ4,λ5
            λ2 = Array{Float64,1}(UndefInitializer(), 5)
            λ2[1] = u2[1]-c2[1]
            λ2[2] = u2[1]
            λ2[3] = u2[1]
            λ2[4] = u2[1]
            λ2[5] = u2[1]+c2[1]

            # get L1 from downstream
            L[1] = λ2[1].*((p2[2]-p2[1])./Δx-ρ2[1].*c2[1].*(u2[2]-u2[1])./Δx)

        return L
        end


"""
    get the characteristic wave amplitudes L from constant enthalpy conditions
"""

    function get_L_from_constant_h(uu::Array,everythinginitial,Δx::Float64)

    gamma=everythinginitial.gamma

    uueverything = UUtoEverything(uu,gamma)

    u = uueverything.u
    ρ = uueverything.ρ
    c = uueverything.c
    p = uueverything.p
    h = uueverything.h

    L = get_L_from_nonreflect(uu,everythinginitial,Δx)
    dhdt = 0
    L[1] = (-dhdt .* ρ[1] .* c[1].^2 ./h[1] + L[2]) .*2 ./(gamma-1) - L[5]

    return L
    end

"""
    get the characteristic wave amplitudes L from constant pressure conditions
"""

    function get_L_from_constant_p(uu::Array,everythinginitial,Δx::Float64)

    L = get_L_from_nonreflect(uu,everythinginitial,Δx)
    dpdt = 0
    L[1] = -L[5] - 2 .* dpdt

    return L
    end

    """
        get the characteristic wave amplitudes L from nonreflect conditions
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
            get the characteristic wave amplitudes L from constant enthalpy inflow
        """

        function get_L_from_inflow_h(uu::Array,everythinginitial,Δx::Float64)

                    # import variables
                    gamma = everythinginitial.gamma


                    uueverything = UUtoEverything(uu,gamma)

                    u = uueverything.u
                    ρ = uueverything.ρ
                    c = uueverything.c
                    p = uueverything.p
                    h = uueverything.h

                    dudt=0
                    dhdt=0

                    # get λ1,λ2,λ3,λ4,λ5
                    λ = Array{Float64,1}(UndefInitializer(), 5)
                    λ[1] = u[end]-c[end]
                    λ[2] = u[end]
                    λ[3] = u[end]
                    λ[4] = u[end]
                    λ[5] = u[end]+c[end]

                    # get L1,L2,L3,L4,L5
                    L = Array{Float64,1}(UndefInitializer(), 5)
                    L[1]=λ[1].*((p[2]-p[1])./Δx-ρ[1].*c[1].*(u[2]-u[1])./Δx)
                    L[5]=L[1] - 2 .* ρ[1] .* c[1] .* dudt
                    L[2]=0.5.*(gamma-1).*(L[5]+L[1]) +ρ[1].*c[1].*c[1]./h[1].*dhdt
                    L[3]=0
                    L[4]=0

                return L
                end

        """
            get the d from L for inflow
        """

        function get_d_from_L_inflow(uu::Array,everythinginitial,L::Array)

                # import variables
                gamma = everythinginitial.gamma

                uueverything = UUtoEverything(uu,gamma)

                u = uueverything.u
                ρ = uueverything.ρ
                c = uueverything.c
                p = uueverything.p
                h = uueverything.h

                dhdt = 0

                # get d1,d2,d3,d4,d5
                d = Array{Float64}(UndefInitializer(), 5)
                d[1] = 1 ./ (c[1].^2).*(L[2]+0.5.*(L[5]+L[1]))
                d[2] = 0.5 .* (L[5]+L[1])
                d[3] = 0.5 ./ρ[1]./c[1] .* (L[5]-L[1])
                d[4] = 0
                d[5] = 0

                return d
                end

        """
            get the d from L for outflow
        """

        function get_d_from_L_outflow(uu::Array,everythinginitial,L::Array)

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
    d = get_d_from_L_outflow(uu,everythinginitial,L)

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

    """
        still working on this
    """

    function set_outlet_costant_p_boundary!(uu::Array,everythinginitial,Δx::Float64,Δt::Float64)

    L = get_L_from_constant_p(uu,everythinginitial,Δx)
    d = get_d_from_L_outflow(uu,everythinginitial,L)

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

    """
        still working on this
    """

    function set_outlet_costant_h_boundary!(uu::Array,everythinginitial,Δx::Float64,Δt::Float64)

    L = get_L_from_constant_h(uu,everythinginitial,Δx)
    d = get_d_from_L_outflow(uu,everythinginitial,L)

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

    """
        still working on this
    """

    function set_inlet_constant_h!(uu::Array,everythinginitial,Δx::Float64,Δt::Float64)

        L = get_L_from_inflow_h(uu,everythinginitial,Δx)
        d = get_d_from_L_inflow(uu,everythinginitial,L)

        gamma = everythinginitial.gamma
        h = everythinginitial.h

        uueverything = UUtoEverything(uu,gamma)

        u = uueverything.u
        ρ = uueverything.ρ

        uubegin=Array{Float64,1}(UndefInitializer(), 3)
        uubegin[1]=uu[1,1] + (-d[1]).*Δt
        uubegin[2]=uu[2,1] + (-u[1].*d[1]-ρ[1].*d[3]+0).*Δt
        uubegin[3]=uu[3,1] + (-0.5 .* u[1].*u[1].*d[1]-d[2]./(gamma-1) - ρ[1].*u[1].*d[3] + 0).*Δt

        return uubegin
        end

        function set_inlet_interface!(uu1::Array,uu2::Array,everythinginitial1,everythinginitial2,Δx::Float64,Δt::Float64)

            L = get_L_from_interface(uu1,uu2,everythinginitial1,everythinginitial2,Δx)
            d = get_d_from_L_inflow(uu2,everythinginitial1,L)

            gamma = everythinginitial2.gamma
            h = everythinginitial2.h

            uueverything = UUtoEverything(uu2,gamma)

            u = uueverything.u
            ρ = uueverything.ρ

            uubegin=Array{Float64,1}(UndefInitializer(), 3)
            uubegin[1]=uu2[1,1] + (-d[1]).*Δt
            uubegin[2]=uu2[2,1] + (-u[1].*d[1]-ρ[1].*d[3]+0).*Δt
            uubegin[3]=uu2[3,1] + (-0.5 .* u[1].*u[1].*d[1]-d[2]./(gamma-1) - ρ[1].*u[1].*d[3] + 0).*Δt

            return uubegin
            end


        function set_outlet_interface!(uu1::Array,uu2::Array,everythinginitial1,everythinginitial2,Δx::Float64,Δt::Float64)

            L = get_L_from_interface(uu1,uu2,everythinginitial1,everythinginitial2,Δx)
            d = get_d_from_L_outflow(uu1,everythinginitial2,L)

            gamma = everythinginitial1.gamma
            h = everythinginitial1.h

            uueverything = UUtoEverything(uu1,gamma)

            u = uueverything.u
            ρ = uueverything.ρ

            uuend=Array{Float64,1}(UndefInitializer(), 3)
            uuend[1]=uu1[1,end] + (-d[1]).*Δt
            uuend[2]=uu1[2,end] + (-u[end].*d[1]-ρ[end].*d[3]+0).*Δt
            uuend[3]=uu1[3,end] + (-0.5 .* u[end].*u[end].*d[1]-d[2]./(gamma-1) - ρ[end].*u[end].*d[3] + 0).*Δt

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
