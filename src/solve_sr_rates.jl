using SpecialFunctions
using Random
using OrdinaryDiffEq
using Statistics
using Distributions
using ForwardDiff: gradient, derivative, Dual, Partials, hessian, jacobian
using DelimitedFiles
using NLsolve
using Interpolations
using HypergeometricFunctions
using SpinWeightedSpheroidalHarmonics
using NPZ
using LinearAlgebra
using WignerSymbols
# NOTE: Core/constants.jl should be included by the parent file
include("heunc.jl")
include("state_utils.jl")

"""
    assoc_legendre_P(l::Int, m::Int, x::Float64) -> Float64

Compute the associated Legendre polynomial P_l^m(x) using recursion.
"""
function assoc_legendre_P(l::Int, m::Int, x::Float64)::Float64
    abs_m = abs(m)

    # Base cases
    if abs_m > l
        return 0.0
    end

    # P_m^m(x) = (-1)^m * (2m-1)!! * (1-x²)^(m/2)
    if l == abs_m
        return (-1)^abs_m * prod(2*k - 1 for k in 1:abs_m) * (1 - x^2)^(abs_m / 2)
    end

    # P_{m+1}^m(x) = x * (2m+1) * P_m^m(x)
    if l == abs_m + 1
        pm_m = assoc_legendre_P(abs_m, abs_m, x)
        return x * (2*abs_m + 1) * pm_m
    end

    # Recurrence: (l-m)*P_l^m(x) = x*(2l-1)*P_{l-1}^m(x) - (l+m-1)*P_{l-2}^m(x)
    p_lm2 = assoc_legendre_P(l - 2, abs_m, x)
    p_lm1 = assoc_legendre_P(l - 1, abs_m, x)
    return (x * (2*l - 1) * p_lm1 - (l + abs_m - 1) * p_lm2) / (l - abs_m)
end

"""
    sphericalY(l::Int, m::Int, theta::Real, phi::Real) -> Complex

Compute the standard spherical harmonic Y_l^m(θ, φ).
Uses closed-form expressions via associated Legendre polynomials.
"""
function sphericalY(l::Int, m::Int, theta::Float64, phi::Float64)::Complex{Float64}
    # Validate quantum numbers
    if abs(m) > l || l < 0
        return 0.0 + 0.0im
    end

    # Normalization constant: sqrt((2l+1)/(4π) * (l-|m|)! / (l+|m|)!)
    # Use BigFloat to avoid factorial overflow for large m
    l_minus_m = l - abs(m)
    l_plus_m = l + abs(m)
    norm_squared = (2*l + 1) / (4*pi) * factorial(big(l_minus_m)) / factorial(big(l_plus_m))
    norm = sqrt(Float64(norm_squared))

    # Associated Legendre polynomial P_l^|m|(cos(θ))
    cos_theta = cos(theta)
    legendre = assoc_legendre_P(l, abs(m), cos_theta)

    # Phase factor: (-1)^m for m < 0
    phase = ifelse(m >= 0, 1.0, (-1.0)^abs(m))

    # Full expression: norm * P_l^|m|(cos θ) * exp(i m φ)
    return phase * norm * legendre * exp(im * m * phi)
end

function ergL(n, l, m, massB, MBH, a; full=true)
    # Key to next level
    alph = GNew * MBH * massB
    if full
        if (l > 0)
            return massB .* (1.0 .- alph.^2 ./ (2 .* n.^2) .- alph.^4 ./ (8 * n.^4) .+ alph.^4 ./ n^4 .* (2 .* l .- 3 .* n .+ 1) ./ (l .+ 0.5) .+ 2 .* a .* m .* alph.^5 ./ n^3 ./ (l .* (l .+ 0.5) .* (l .+ 1)))
        else
            return massB .* (1.0 .- alph.^2 ./ (2 .* n.^2) .- alph.^4 ./ (8 * n.^4) .+ alph.^4 ./ n^4 .* (2 .* l .- 3 .* n .+ 1) ./ (l .+ 0.5))
        end
    else
        return massB .* (1.0 .- alph.^2 ./ (2 .* n.^2) .- alph.^4 ./ (8 * n.^4))
    end
end


function sr_rates(n, l, m, massB, MBH, aBH; impose_low_cut=0.001, solve_322=true)

    alph = GNew .* MBH .* massB

    rP = nothing
    if aBH .> maxSpin
        rP = 1.0 .+ sqrt.(1 - maxSpin .^2)
        aBH = maxSpin
    elseif aBH .< 0.0
        rP = 2.0
        aBH = 0.0
    else
        rP = 1.0 .+ sqrt.(1 - aBH.^2)
    end
    rP *= (GNew .* MBH)
    OmegaH = aBH ./ (2 .* (GNew .* MBH) .* (1 .+ sqrt.(1 .- aBH.^2)))
    Anl = 2 .^(4 .* l .+ 1) .* factorial(big(Int(l .+ n))) ./ (n.^(2 .* l .+ 4) .* factorial(big(n .- l .- 1)))
    Anl *= (factorial(big(Int(l))) ./ (factorial(big(Int(2 .* l))) .* factorial(big(Int(2 .* l .+ 1))))).^2
    Chilm = 1.0
    for k in 1:Int(l)
        Chilm *= (k.^2 .* (1.0 .- aBH.^2) .+ (aBH * m .- 2 .* (rP ./ (GNew .* MBH)) .* alph).^2)
    end
    Gamma_nlm = 2 * massB .* rP .* (m .* OmegaH .- ergL(n, l, m, massB, MBH, aBH)) .* alph.^(4 .* l + 4) .* Anl .* Chilm
    

    if Gamma_nlm > 0.0
        return Gamma_nlm
    elseif (Gamma_nlm < 0.0)&&(m .* OmegaH .> ergL(n, l, m, massB, MBH, aBH))
        return 0.0
    else
        return Gamma_nlm
    end
end

function radial_bound_NR(n, l, m, mu, M, r; physU=false)
    # r in units r -> r / (GM)
    alph = GNew * M * mu
    
    if !physU
        a0 = 1 / alph.^2
    else
        a0 = 1 / (mu * alph)
    end
    rF = sqrt.((2 / (n * a0)).^3 * factorial(big(n - l - 1)) / (2 * n * factorial(big(n + l)))) .* exp.( - r ./ (n .* a0)) .* (2 .* r ./ (n .* a0)).^l .* generalized_laguerre(n - l - 1, 2 * l + 1, 2 .* r ./ (n .* a0)) # I've thrown through a factor of (GNew * M)^{3/2}
    return rF
end

function radial_inf_NR(k, l, mu, M, r)
    # r assumed to be in [G M]
    # k unitless
    alph = GNew * M * mu
    a0 = 1 / alph.^2
    rF = zeros(Complex, length(r))
    for i in 1:length(r)
        rF[i] = 2 .* k .* exp.(pi ./ (2 .* k .* a0)) .* abs.(gamma(l + 1 - im ./ (k .* a0))) ./ factorial(big(2 * l + 1)) .* (2 .* k .* r[i]).^l .* exp.(-im * k * r[i]) .* pFq((im ./ (k .* a0) + l + 1,), (2 * l + 2, ), 2 * im .* k .* r[i])
    end
    return rF # unitless
end


function freq_shifts(mu, M, a, n1, l1, m1, n2, l2, m2;  rpts=500, rmaxT=100, Nang=100000, epsil_2=1.0, Npts_Bnd=500)
    rp = 1 + sqrt.(1 - a.^2)
    alph = mu * GNew * M
    

    erg_1G = ergL(n1, l1, m1, mu, M, a)
    q = abs.(real(- sqrt.(alph.^2 .- (erg_1G .* GNew .* M).^2)))
    rmax = rmaxT ./ q
    rlist = 10 .^(range(log10.(rp), log10.(rmax), rpts))
    
    rl, r1, erg_1 = solve_radial(mu, M, a, n1, l1, m1; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true)
    itp = LinearInterpolation(log10.(rl), log10.(r1), extrapolation_bc=Line())
    rf_1 = 10 .^itp(log10.(rlist))
    Z1 = spheroidals(l1, m1, a, erg_1 ./ (GNew .* M))
    
    # rf_1 = radial_bound_NR(n1, l1, m1, mu, M, rlist)

    
    if (n2 == n1)&&(l2==l1)&&(m2==m1)
        rf_2 = rf_1
        erg_2 = erg_1
        Z2 = Z1
        mult_fac = 1.0 ./ 2.0
    else
        rl, r2, erg_2 = solve_radial(mu, M, a, n2, l2, m2; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true)
        itp = LinearInterpolation(log10.(rl), log10.(r2), extrapolation_bc=Line())
        rf_2 = 10 .^itp(log10.(rlist))
        Z2 = spheroidals(l2, m2, a, erg_2 ./ (GNew .* M))
        mult_fac = 1.0
    end
    
    function func_ang(x)
        return real(Z1.(x[1], x[2]) .* Z2.(x[1], x[2]) .* conj(Z1.(x[1], x[2])) .* conj(Z2.(x[1], x[2])))
    end
    
    
    thetaV = acos.(1.0 .- 2.0 .* rand(Nang))
    phiV =  2.0 * pi .* rand(Nang)
    CG = 0.0
    for i in 1:Nang
        CG += func_ang([thetaV[i], phiV[i]])
    end
    CG *= 4 * pi / Nang
    
    r_integrd = rf_1 .* rf_2 .* conj(rf_1) .* conj(rf_2) .* rlist.^2
    radial_int = trapz(r_integrd, rlist) ./ (GNew * M).^3
    
    delt_om = -1.0 ./ (4 * mu^2) .* radial_int .* CG .* mult_fac
    
    return Float64(real(alph.^2 .* delt_om .* epsil_2)) # \delta Omega [Needs to be multiplied by (M_pl / f_a)^2], do dF / (alpha^5 * mu) to get coefficient hitting epsilon_x

end

function s_rate_bnd(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; kpts=10, rpts=2000, rmaxT=100, inf_nr=true, Nang=100000, Npts_Bnd=1000, debug=false, include_cont=true, Ntot_safe=5000, sve_for_test=false, bnd_thresh=1e-3, use_analytic=false, eps_r=1e-10, NON_REL=true)
    
    rp = 1 + sqrt.(1 - a.^2)
    alph = mu * GNew * M
    
    # approx energy estimates
    erg_1G = ergL(n1, l1, m1, mu, M, a)
    erg_2G = ergL(n2, l2, m2, mu, M, a)
    erg_3G = ergL(n3, l3, m3, mu, M, a)
    
    q1 = abs.(real(- sqrt.(alph.^2 .- (erg_1G .* GNew .* M).^2)))
    q2 = abs.(real(- sqrt.(alph.^2 .- (erg_2G .* GNew .* M).^2)))
    q3 = abs.(real(- sqrt.(alph.^2 .- (erg_3G .* GNew .* M).^2)))
    
    qmax = maximum([q1, q2, q3])
    qmin = minimum([q1, q2, q3])
    rmax = rmaxT ./ qmax
    rlist = 10 .^(range(log10.(rp .* (1.0 .+ eps_r)), log10.(rmax), rpts))
    
    
 
    if !use_analytic || !NON_REL
        rl, r1, erg_1 = solve_radial(mu, M, a, n1, l1, m1; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
        itp = LinearInterpolation(log10.(rl), r1, extrapolation_bc=Line())
        rf_1 = itp(log10.(rlist))
        # rf_1 = radial_bound_NR(n1, l1, m1, mu, M, rlist)
        # erg_1 = erg_1G .* GNew * M
    else
        rf_1 = radial_bound_NR(n1, l1, m1, mu, M, rlist)
        erg_1 = erg_1G .* GNew * M
    end
    

    
    if (n2 == n1)&&(l2==l1)&&(m2==m1)
        rf_2 = rf_1
        erg_2 = erg_1
    else
        if !use_analytic || !NON_REL
            rl, r2, erg_2 = solve_radial(mu, M, a, n2, l2, m2; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
            itp = LinearInterpolation(log10.(rl), r2, extrapolation_bc=Line())
            rf_2 = itp(log10.(rlist))
            # rf_2 = radial_bound_NR(n2, l2, m2, mu, M, rlist)
            # erg_2 = erg_2G .* GNew * M
        else
            rf_2 = radial_bound_NR(n2, l2, m2, mu, M, rlist)
            erg_2 = erg_2G .* GNew * M
        end
    end

    if !use_analytic || !NON_REL
        rl, r3, erg_3 = solve_radial(mu, M, a, n3, l3, m3; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
        itp = LinearInterpolation(log10.(rl), r3, extrapolation_bc=Line())
        rf_3 = itp(log10.(rlist))
        # rf_3 = radial_bound_NR(n3, l3, m3, mu, M, rlist)
        # erg_3 = erg_3G .* GNew * M
    else
        rf_3 = radial_bound_NR(n3, l3, m3, mu, M, rlist)
        erg_3 = erg_3G .* GNew * M
    end
    

    erg_ind = erg_1 .+ erg_2 - erg_3
    k_ind_2 = mu.^2 .- (erg_ind ./ (GNew .* M)).^2
    

    Z1 = spheroidals(l1, m1, a, erg_1 ./ (GNew .* M))
    Z2 = spheroidals(l2, m2, a, erg_2 ./ (GNew .* M))
    Z3 = spheroidals(l3, m3, a, erg_3 ./ (GNew .* M))
    
    
        

    # compute bound contribution
    kmin = 0.01 .* alph^2
    kmax = 2.0 .* alph^2
    kk_list = 10 .^LinRange(log10.(kmin), log10.(kmax),  kpts)
   
    ck_list = zeros(Complex, kpts)
 
    UnF = GNew .* M
    bound_c = 0.0
    contin_c = 0.0
    done_nmax = false
    n = 1
    
    function func_ang(x, Zf1, Zf2, Zf3, Zf4)
        return real(Zf1.(x[1], x[2]) .* Zf2.(x[1], x[2]) .* conj(Zf3.(x[1], x[2])) .* conj(Zf4.(x[1], x[2])))
    end
        
    if !NON_REL
        rmm = 1.0 - sqrt.(1.0 - a.^2)
        rout_star = rlist .+ 2.0 .* rp ./ (rp .- rmm) .* log.(rlist ./ rp .- 1.0) .- 2.0 .* rmm .* log.(rlist ./ rmm .- 1.0) ./ (rp .- rmm)
        itp_rrstar = LinearInterpolation(rout_star, rlist, extrapolation_bc=Line())
        itp_rrstar_inv = LinearInterpolation(rlist, rout_star, extrapolation_bc=Line())
    end
        
    while !done_nmax
        
        rl, r4, erg_4 = solve_radial(mu, M, a, n, 0, 0; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
        if NON_REL
            r4 = radial_bound_NR(n, 0, 0, mu, M, rlist)
            rl = rlist
            erg_4G = alph .* (1 - alph.^2 ./ (2 .* n.^2))
            # Use tolerance-based comparison instead of exact equality for floating-point numbers
            if abs(erg_4G - erg_ind) > 1e-10 * abs(erg_ind)
                erg_4 = erg_4G
            end
        end
        if NON_REL
            itp = LinearInterpolation(log10.(rl), r4, extrapolation_bc=Line())
            rf_4 = itp(log10.(rlist))
 
        else
            itp = LinearInterpolation(log10.(rl), (r4 .* sqrt.(rl.^2 .+ a.^2)), extrapolation_bc=Line())
            rf_4 = itp(log10.(rlist))
        end
        
        
        Z4 = spheroidals(0, 0, a, erg_4 ./ (GNew .* M))
        
        thetaV = acos.(1.0 .- 2.0 .* rand(Nang))
        phiV = rand(Nang) .* 2*pi
        CG = 0.0
        for i in 1:Nang
            CG += func_ang([thetaV[i], phiV[i]], Z1, Z2, Z3, Z4)
        end
        CG *= 4*pi / Nang
       
        if NON_REL
            r_integrd = (rf_1 ./ UnF^(3/2)) .* (rf_2 / UnF^(3/2)) .* conj(rf_3 / UnF^(3/2)) .* conj(rf_4 / UnF^(3/2)) .* (rlist * UnF).^2
            radial_int = trapz(r_integrd, rlist * UnF)
            
            kdiff_sq = (erg_4.^2 - erg_ind^2) ./ (GNew .* M).^2
        
            ff = CG .* radial_int ./ (2 .* mu).^(3 / 2) ./ (kdiff_sq)
            ### double factor 1/2
            if (n1==n2)&&(l1==l2)&&(m1==m2)
                ff /= 2
            end
            res_n = (rf_4[1] ./ UnF^(3/2)) .* ff
        else
            # dropping a^2 term in CG_2
            r_integrd = rf_1 .* rf_2 .* conj(rf_3) .* conj(rf_4) .* (rlist.^2 .* (rlist.^2 .- 2 .* rlist .+ a.^2) ./ (rlist.^2 .+ a.^2).^(3/2)) .* rlist.^2 ./ (rlist.^2 .+ a.^2)
            radial_int = trapz(r_integrd, rlist)
            
            kdiff_sq = (erg_4.^2 - erg_ind^2)
        
            ff = CG .* radial_int ./ (2 .* alph).^(3 / 2) ./ (kdiff_sq)
            ### double factor 1/2
            if (n1==n2)&&(l1==l2)&&(m1==m2)
                ff /= 2
            end
            res_n = (rf_4[1] ./ sqrt.(rlist[1].^2 .+ a.^2)) .* ff ./ (GNew .* M)
        end
        

        if debug
            print(n, "\t", Float64(abs.(res_n)), "\n")
        end
        
        bound_c += res_n
        if abs.(res_n / bound_c) < bnd_thresh
            done_nmax = true
        end
        n += 1
        
        done_nmax = true ### TEMP HOLD
    end
    


    if include_cont && NON_REL
        # compute continuous contribution
        for i in 1:kpts
            k = kk_list[i] ./ (GNew .* M) # physical units
            erg_New = sqrt.(k.^2 .+ mu.^2)
            Z4 = spheroidals(0, 0, a, erg_New)
        
            thetaV = acos.(1.0 .- 2.0 .* rand(Nang))
            phiV = rand(Nang) .* 2*pi
            CG = 0.0
            for i in 1:Nang
                CG += func_ang([thetaV[i], phiV[i]], Z1, Z2, Z3, Z4)
            end
            CG *= 4*pi / Nang
        
            if inf_nr
                out_goingR = radial_inf_NR(k .* GNew .* M, 0, mu, M, rlist)
            else
                rl, r4 = radial_inf(erg_New .* GNew .* M, mu, M, a, 0, 0; rmax_val=rmax, rpts=rpts)
                # itp = LinearInterpolation(log10.(rl[3:end]), log10.(r4[3:end]), extrapolation_bc=Line())
                if NON_REL
                    itp = LinearInterpolation(log10.(rl[3:end]), r4[3:end], extrapolation_bc=Line())
                else
                    itp = LinearInterpolation(log10.(rl[3:end]), (r4[3:end] .* sqrt.(rl[3:end].^2 .+ a.^2)), extrapolation_bc=Line())
                end
                out_goingR = itp(log10.(rlist))
                
                
            end
            
            if NON_REL
                r_integrd = (rf_1 / UnF^(3/2)) .* (rf_2 / UnF^(3/2)) .* conj(rf_3 / UnF^(3/2)) .* conj(out_goingR ./ UnF) .* (rlist * UnF).^2
                radial_int = trapz(r_integrd, rlist * UnF)
            else
                r_integrd = rf_1 .* rf_2 .* conj(rf_3) .* conj(out_goingR .* sqrt.(rlist.^2 .+ a.^2)) .* (rlist.^2 .* (rlist.^2 .- 2 .* rlist .+ a.^2) ./ (rlist.^2 .+ a.^2).^(3/2)) .* rlist.^2 ./ (rlist.^2 .+ a.^2)
                radial_int = trapz(r_integrd, rlist) ./ UnF.^(5 / 2)
            end
        
       
            ff = CG .* radial_int ./ (2 .* mu).^(3 / 2) ./ (k_ind_2 .+ k.^2) ./ (2 * pi)
        
            ### double factor 1/2
            if (n1==n2)&&(l1==l2)&&(m1==m2)
                ff /= 2
            end
  
            
            if !isnan.(out_goingR[1] .* ff ./ UnF)
                ck_list[i] = out_goingR[1] .* ff ./ UnF
            else
                print("getting nan... \t", i, "\t", out_goingR[1:3], "\t", ff, "\n")
            end
            if debug
                print("cont \t", i, "\t", kk_list[i], "\t", Float64(abs.(ck_list[i])), "\t", ff, "\n")
            end
        end
        contin_c += trapz(ck_list, kk_list / UnF)
    end
    
    psi_1 = bound_c .+ contin_c
    print("bound contribution:  ", Float64(abs.(bound_c)), "\t contin: ", Float64(abs.(contin_c)), "\n")
        
    
    lam = (mu ./ (M_pl .* 1e9))^2
    rate_out = 4 .* alph.^2 .* (1 .+ sqrt.(1 - a.^2)) .* Float64(real(psi_1 .* conj(psi_1))) .* lam^2
    
    return rate_out ./ mu^2 .* (GNew * M^2 * M_to_eV)^2 # unitless [gamma / mu]
end

function s_rate_inf(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3, lF_min; rpts=4000, rmaxT=90,  sve_for_test=false, inf_nr=false, Npts_Bnd=1000, Nang=300000, debug=false, Ntot_safe=2000, xtol=1e-2, ftol=1e-2, iter=20)
    # lF_min is l relative to minimal value needed for non-zero m
    # returns scattering re-normalized \gamma (basically just radial integral ratio)
    

    
    rp = (1.0 .+ sqrt.(1 - a.^2))
    alph = mu .* GNew .* M
    
    # rough estimate energies
    erg_1 = ergL(n1, l1, m1, mu, M, a)
    erg_2 = ergL(n2, l2, m2, mu, M, a)
    erg_3 = ergL(n3, l3, m3, mu, M, a)
    erg_New = erg_1 .+ erg_2 - erg_3
    # k = sqrt.(erg_New.^2 .- mu.^2)
    if debug
        print("erg new \t", erg_New .* GNew * M, "\n")
    end
    k = sqrt.((erg_New .* GNew * M).^2 .- alph.^2)
    mF = (m1 .+ m2 - m3)
    lF = abs.(mF) + lF_min
    
    
    
    q1 = abs.(real(- sqrt.(alph.^2 .- (GNew .* M .* erg_1).^2)))
    q2 = abs.(real(- sqrt.(alph.^2 .- (GNew .* M .* erg_2).^2)))
    q3 = abs.(real(- sqrt.(alph.^2 .- (GNew .* M .* erg_3).^2)))
    q = maximum([q1, q2, q3])
    qcheck = 6 ./ (k .* GNew .* M)
    rmax = rmaxT ./ q
    
    rlist = 10 .^(range(log10.(rp), log10.(rmax), rpts))
    # print("New erg \t", erg_New .* GNew .* M, "\n")


    Z1 = spheroidals(l1, m1, a, erg_1)
    Z2 = spheroidals(l2, m2, a, erg_2)
    Z3 = spheroidals(l3, m3, a, erg_3)
    Z4 = spheroidals(lF, mF, a, erg_New)
   
    thetaV = acos.(1.0 .- 2.0 .* rand(Nang))
    phiV = rand(Nang) .* 2*pi

    function func_ang(x)
        return real(Z1.(x[1], x[2]) .* Z2.(x[1], x[2]) .* conj(Z3.(x[1], x[2])) .* conj(Z4.(x[1], x[2]))) #
    end
    CG = 0.0
    ang_fac2 = 0.0
    
    for i in 1:Nang
        CG += func_ang([thetaV[i], phiV[i]])
    end
    CG *= 4*pi / Nang
   
    
    rl, r1 = solve_radial(mu, M, a, n1, l1, m1; rpts=Npts_Bnd, rmaxT=rmaxT, Ntot_safe=Ntot_safe)
    itp = LinearInterpolation(log10.(rl), r1, extrapolation_bc=Line())
    rf_1 = itp(log10.(rlist))
    
    if sve_for_test
        writedlm("test_store/check_real_0.dat", hcat(rl, float(real(r1))))
        writedlm("test_store/check_imag_0.dat", hcat(rl, float(imag(r1))))
    end
    
    if (n2 == n1)&&(l2==l1)&&(m2==m1)
        rf_2 = rf_1
    else
        rl, r2 = solve_radial(mu, M, a, n2, l2, m2; rpts=Npts_Bnd, rmaxT=rmaxT, Ntot_safe=Ntot_safe)
        itp = LinearInterpolation(log10.(rl), r2, extrapolation_bc=Line())
        rf_2 = itp(log10.(rlist))
    end

    
    rl, r3 = solve_radial(mu, M, a, n3, l3, m3; rpts=Npts_Bnd, rmaxT=rmaxT, Ntot_safe=Ntot_safe)
    itp = LinearInterpolation(log10.(rl), r3, extrapolation_bc=Line())
    rf_3 = itp(log10.(rlist))
    
    
    if inf_nr
        out_goingR = radial_inf_NR(k, lF, mu, M, rlist)
    else
        rl, r4 = radial_inf(erg_New .* GNew .* M, mu, M, a, lF, mF; rmax_val=rmax, rpts=rpts, debug=debug, xtol=xtol, ftol=ftol, iter=iter)
        # itp = LinearInterpolation(log10.(rl[3:end]), (r4[3:end]), extrapolation_bc=Line())
        # out_goingR = itp(log10.(rlist))
        itp = LinearInterpolation(log10.(rl[3:end]), r4[3:end], extrapolation_bc=Line())
        out_goingR = itp(log10.(rlist))
    end
    
    if sve_for_test
        writedlm("test_store/check_real_1.dat", hcat(rlist, float(real(rf_1))))
        writedlm("test_store/check_imag_1.dat", hcat(rlist, float(imag(rf_1))))
        writedlm("test_store/check_real_2.dat", hcat(rlist, float(real(rf_2))))
        writedlm("test_store/check_imag_2.dat", hcat(rlist, float(imag(rf_2))))
        writedlm("test_store/check_real_3.dat", hcat(rlist, float(real(rf_3))))
        writedlm("test_store/check_imag_3.dat", hcat(rlist, float(imag(rf_3))))
        writedlm("test_store/check_real_4.dat", hcat(rlist, float(real(out_goingR))))
        writedlm("test_store/check_imag_4.dat", hcat(rlist, float(imag(out_goingR))))
    end
    
    unitsF = GNew * M
    r_integrd = (rf_1 ./ unitsF^(3/2))  .* (rf_2 ./ unitsF^(3/2)) .* conj(rf_3 / unitsF^(3/2)) .* conj(out_goingR / unitsF) .* (rlist * unitsF).^2
    radial_int = trapz(r_integrd, rlist * unitsF)
    
    # Rbound [1/GM]^{3/2}, R_c [eV], r [GM]
    ff = CG .* 4 * pi * (- im).^lF ./ (2 .* k ./ unitsF) .* radial_int ./ (2 .* mu).^(3 / 2)
    
    ### double factor 1/2
    if (n1==n2)&&(l1==l2)&&(m1==m2)
        ff /= 2
    end
    
    lam = (mu ./ (M_pl .* 1e9))^2
    rate_out = Float64(real(2 .* erg_New * (k / unitsF) ./ (4 * pi)^2 .* ff .* conj(ff) .* lam.^2))
    
    return rate_out ./ mu^2 .* (GNew * M^2 * M_to_eV)^2 # unitless [gamma / mu]
    # Return: rate  (evaluated at f_a = Mpl)

    
end

function solve_radial(mu, M, a, n, l, m; rpts=1000, rmaxT=50, debug=false, iter=500, xtol=1e-50, ftol=1e-90, sve=false, fnm="test_store/WF_", return_erg=false, Ntot_safe=5000, eps_r=1e-10, QNM=false, QNM_ergs=nothing, pre_compute_erg=nothing, prec=300, use_heunc=true)
    ### dolan 2007
    
    
    setprecision(BigFloat, prec)
    
    alph = GNew .* M .* mu
    
    if m == 1
        Ntot = 1500
    elseif (m==2)
        Ntot = 4250
    elseif (m==3)
        Ntot = 5000
    else
        Ntot = Ntot_safe
    end
    
    if a > 0.95
        Ntot *= 2.5
        Ntot = Int(Ntot)
    end
    if Ntot_safe > Ntot
        Ntot = Ntot_safe
    end
    if !QNM
        if isnothing(pre_compute_erg)
            wR, wI = find_im_part(mu, M, a, n, l, m; debug=debug, iter=iter, xtol=xtol, ftol=ftol, return_both=true, for_s_rates=true, QNM=false, Ntot_force=Ntot_safe)
        else
            wR = BigFloat(real(pre_compute_erg))
            wI = BigFloat(imag(pre_compute_erg))
        end
    else
        wR = real(QNM_ergs)
        wI = imag(QNM_ergs)
    end
    erg = wR .+ im .* wI
    if debug
        print("erg \t", erg, "\n")
    end
    
    
    q = - sqrt.(alph.^2 .- erg.^2)
    if !QNM
        if real(q) > 0
            q *= -1
        end
    else
        if real(q) < 0
            q *= -1
        end
    end
    

    # rmax = rmaxT * 1.0 ./ abs.(real(q))
    rmax = Float64.(100 * n ./ alph.^2)
    # print("Check \t", rmax, "\t", rmaxT * 1.0 ./ abs.(real(q)), "\n")
    
    if !use_heunc
        rm = (1.0 .- sqrt.(1 - a.^2))
        rp = (1.0 .+ sqrt.(1 - a.^2))
        OmH = a  ./ (2 .* rp)

        
        sigm = 2 .* rp .* (erg .- OmH .* m) ./ (rp .- rm)
        chi = (alph.^2 .- 2 .* erg.^2) ./ q
        b = sqrt.(1 - a.^2)
        
        gam = im * a * sqrt.(erg.^2 .- alph.^2)
        cc2 = a.^2 .* (erg.^2 .- alph.^2)
        LLM = l * (l + 1)
        # LLM += 2 .* cc2 .* (m.^2 .- l .* (l + 1) + 0.5) ./ ((2 * l - 1) .* (2 * l + 3)) # alternative describ 2nd order
        LLM += (-1 + 2 * l * (l + 1) - 2 * m.^2) * gam.^2 ./ (-3 + 4 * l * (l + 1))
        LLM += ((l - m - 1 * (l - m) * (l + m) * (l + m - 1)) ./ ((-3 + 2 * l) * (2 * l - 1).^2) - (l + 1 - m) * (2 * l - m) * (l + m + 1) * (2 + l + m) ./ ((3 + 2 * l).^2 * (5 + 2 * l))) * gam.^4 ./ (2 * (1 + 2 * l))
        LLM += (4 * ((-1 + 4 * m^2) * (l * (1 + l) * (121 + l * (1 + l) * (213 + 8 * l * (1 + l) * (-37 + 10 * l * (1 + l)))) - 2 * l * (1 + l) * (-137 + 56 * l * (1 + l) * (3 + 2 * l * (1 + l))) * m^2 + (705 + 8 * l * (1 + l) * (125 + 18 * l * (1 + l))) * m^4 - 15 * (1 + 46 * m^2))) * gam^6) / ((-5 + 2 * l) * (-3 + 2 * l) * (5 + 2 * l) * (7 + 2 * l) * (-3 + 4 * l * (1 + l))^5)
        
        c0 = 1.0 .- 2.0 * im * erg - 2 * im ./ b .* (erg .- a .* m ./ 2.0)
        c1 = -4.0 .+ 4 * im * (erg - im * q * (1.0 + b)) + 4 * im / b * (erg - a * m / 2.0) .- 2.0 * (erg.^2 .+ q.^2) ./ q
        c2 = 3.0 - 2 * im * erg - 2.0 * (q.^2 - erg.^2) ./ q - 2.0 * im / b * (erg - a * m ./ 2)
        c3 = 2.0 * im * (erg - im * q).^3 ./ q .+ 2 * (erg .- im * q).^2 .* b + q.^2 * a.^2 .+ 2 * im * q * a * m - LLM - 1 - (erg - im * q).^2 ./ q .+ 2 * q * b + 2 * im / b * ( (erg - im * q).^2 ./ q + 1.0) * (erg .- a * m / 2)
        c4 = (erg .- im * q).^4 ./ q.^2 .+ 2 * im * erg * (erg .- im * q).^2 ./ q .- 2 * im ./ b .* (erg .- im * q).^2 ./ q .* (erg .- a * m ./ 2)
        
        function alphaN(nn)
            return nn.^2 .+ (c0 .+ 1) * nn + c0
        end
        function betaN(nn)
            return -2 * nn.^2 .+ (c1 .+ 2) * nn + c3
        end
        function gammaN(nn)
            return nn.^2 .+ (c2 .- 3) * nn + c4
        end
        
        
        rlist = 10 .^(range(log10.(rp  .* (1.0 .+ eps_r)), log10.(rmax), rpts))
        if (wR == 0)&&(wI == 0)
            if return_erg
                return rlist, zeros(length(rlist)), wR .+ im .* wI
            else
                return rlist, zeros(length(rlist))
            end
        end
        
        preF = (rlist .- rp).^(- im .* sigm) .* (rlist .- rm).^(im .* sigm .+ chi .- 1.0) .* exp.(q .* rlist)

        
        nfactor = ((rlist .- rp) ./ (rlist .- rm))
        Rout = zeros(Complex, length(rlist))
        bk = zeros(Complex, Ntot)
        bk[1] = 1.0
        bk[2] = -betaN(0) ./ alphaN(0) .* bk[1]
        
        for i in 2:(Ntot-1)
            bk[i+1] = (-betaN(i-1) .* bk[i] .- gammaN(i-1) .* bk[i-1]) ./ alphaN(i-1)
        end
        
        for i in 1:Ntot
            Rout .+= preF .* bk[i] .* nfactor.^(i-1)
        end
        
       
        
        if QNM
            # nm = 1.0
            alp_trans = 2 .* rp ./ (rp .- rm) .* log.(abs.(rlist .- rp)) .- 2 .* rm ./ (rp .- rm) .* log.(abs.(rlist .- rm))
            Rout .*= exp.(im * erg .* alp_trans) # change coord
            # Rout .*= exp.(-im * erg .* alp_trans) # change coord
        end
     
        nm = trapz(Rout .* conj(Rout) .* rlist.^2, rlist)
        if QNM
            nm = 1.0
        end
        if sve
            normOut = Rout .* conj(Rout) ./ nm
            writedlm(fnm*"_$(alph)_n_$(n)_l_$(l)_m_$(m)_.dat", hcat(float(real(rlist)), float(real(normOut))))
        else
            if return_erg
                return rlist, Rout ./ sqrt.(nm), erg # erg normalized by GM
            else
                return rlist, Rout ./ sqrt.(nm)
            end
        end
    else
        rminus = M - sqrt(M^2 - a^2)
        rplus = M + sqrt(M^2 - a^2)
        function R(r, ν_val)
            return ((2*r*M*alph^2)/(l0 + n0 + 1))^l0 *
                    exp(-((r*M*alph^2)/(l0 + n0 + 1))) *
                    generalized_laguerre(n - l0 - 1, 2 * l0 + 1, (2*r*M*alph^2)/(l0 + n0 + 1))
        end
        function calc_ω_inv(ω)
            return sqrt.(((1 - (ω ./ alph).^2 ) ./ (alph^2*M^2) ).^-1)
        end
        function calc_Pplus(ω_val)
            return (m*a - 2*ω_val*M*rplus)/(rplus - rminus)
        end
        function calc_Pminus(ω_val)
            return (m*a - 2*ω_val*M*rminus)/(rplus - rminus)
        end
        function calc_Aplus(ω_val)
            return 2*((m*a - 2*ω_val*M^2)/(rplus - rminus))^2 -
                   4*ω_val^2*M^2 + (alph^2 - 2*ω_val^2)*M*(rplus - rminus) +
                   (2*M^2 - a^2)*(alph^2 - ω_val^2)
        end
        function calc_Aminus(ω_val)
            return 2*((m*a - 2*ω_val*M^2)/(rplus - rminus))^2 -
                   4*ω_val^2*M^2 - (alph^2 - 2*ω_val^2)*M*(rplus - rminus) +
                   (2*M^2 - a^2)*(alph^2 - ω_val^2)
        end

        ω = erg
        nuV = calc_ω_inv(ω)
        kk = sqrt.(ω.^2 .- alph.^2)
        
        l0 = l
        n0 = n - l0 - 1 # field overtone number
        # Calculate Pplus for the current ν value
        Pplus = calc_Pplus(ω)
        Pmns = calc_Pminus(ω)
        beta = 2 * im .* Pplus
        gam = -2 .* im .* Pmns
        alpha_other = 2 * im * kk .* (rplus .- rminus)
        deltt = calc_Aplus(ω) .- calc_Aminus(ω)
        gam_llm = im * a * sqrt.(ω.^2 .- alph.^2)
        LLM = l0 * (l0 + 1)
        LLM += (-1 + 2 * l0 * (l0 + 1) - 2 * m.^2) * gam_llm.^2 ./ (-3 + 4 * l0 * (l0 + 1))
        LLM += ((l0 - m - 1 * (l0 - m) * (l0 + m) * (l0 + m - 1)) ./ ((-3 + 2 * l0) * (2 * l0 - 1).^2) - (l0 + 1 - m) * (2 * l0 - m) * (l0 + m + 1) * (2 + l0 + m) ./ ((3 + 2 * l0).^2 * (5 + 2 * l0))) * gam_llm.^4 ./ (2 * (1 + 2 * l0))
        LLM += (4 * ((-1 + 4 * m^2) * (l0 * (1 + l0) * (121 + l0 * (1 + l0) * (213 + 8 * l0 * (1 + l0) * (-37 + 10 * l0 * (1 + l0)))) - 2 * l0 * (1 + l0) * (-137 + 56 * l0 * (1 + l0) * (3 + 2 * l0 * (1 + l0))) * m^2 + (705 + 8 * l0 * (1 + l0) * (125 + 18 * l0 * (1 + l0))) * m^4 - 15 * (1 + 46 * m^2))) * gam_llm^6) / ((-5 + 2 * l0) * (-3 + 2 * l0) * (5 + 2 * l0) * (7 + 2 * l0) * (-3 + 4 * l0 * (1 + l0))^5)
        eta = - (Pplus.^2 .+ Pmns.^2 .+ LLM .+ alph.^2 .* rplus.^2 .- ω.^2 .* (4 .* M.^2 .+ 2 .* M .* rplus .+ rplus.^2))
        phi = (alpha_other .- beta .- gam .+ alpha_other .* beta .- beta .* gam) ./ 2 .- eta
        nuN = (alpha_other .+ beta .+ gam .+ alpha_other .* gam  .+ beta .* gam) ./ 2 .+ eta .+ deltt
        
        ### choose r values i want...
        r_max_shrt = 2.0 .^(2.0 .* n .- 2 .* (1 .+ n)) .* gamma(2 .+ 2 .* n) ./ alph.^2 ./ factorial(big(2 .* n - 1)) .* 3.0
        r_max = n ./ alph.^2 .* 100.0
        r_vals = rlist = 10 .^(range(log10.(rplus  .* (1.0 .+ 1e-3)), log10.(r_max), rpts))
        
        rout_temp = R.(r_vals, nuV)
        rout_temp = abs.(rout_temp)
        rout_temp ./= trapz(rout_temp.^2 .* r_vals.^2, r_vals)
        r_thresh = nothing
        for i in 1:length(r_vals)
            if r_vals[end - i + 1] .< r_max_shrt
                r_thresh = r_vals[end - i + 1]
                break
            end
        end
        zz = -(r_vals .- rplus) ./ (rplus .- rminus)

        phinu = Float64(real.(phi + nuN)) .+ im .* Float64(imag.(phi + nuN))
        
        zz_short = zz[abs.(zz) .< r_thresh] # cant call large zz values of heun (takes forever...)
        heunc = solve_heun_switch(phi, phi + nuN, 1.0 + beta, 1.0 + gam, alpha_other, zz_short)
        y_values_short = exp.(-im .* kk .* (r_vals[abs.(zz) .< r_thresh] .- rplus)) .* zz_short.^(im .* Pplus) .* (zz_short .- 1.0).^(- im .* Pmns) .* heunc

        
        y_values = []
        rmax_cut = r_vals[abs.(zz) .< r_thresh][end]
        
        # work backward to find where the solution fails
        imax = length(zz_short)
        if abs.(y_values_short[end]) .> abs.(y_values_short[end-1])
            for i in 1:(length(zz_short)-1)
                if abs.(y_values_short[end-i+1]) .> abs.(y_values_short[end-i])
                    y_values_short[end-i+1] = 0.0 + im * 0.0
                else
                    imax = length(zz_short) - i - 1
                    break
                end
            end
        end
        # now go a little further back because you didn't go far enough
        imax = argmin(abs.(r_vals .- r_vals[imax] .* 0.9))
        r_max_new = r_vals[imax]
        
        for i in 1:rpts
            if r_vals[i] <= r_max_new
                push!(y_values, y_values_short[i])
            else
                push!(y_values, y_values_short[imax] .* exp.(im .* kk .* (r_vals[i] .- r_max_new)) .* (r_vals[i] ./ r_max_new).^2 )
            end
        end
        rlist = r_vals ./ M
        
        
        nm2 = trapz(y_values .* conj.(y_values) .* rlist.^2, rlist)
        y_values ./= sqrt.(nm2)
        
        number_max_cnt = count_local_maxima(real.(y_values .* conj.(y_values)))
        if number_max_cnt != (n - l0)
            println("Nmax count :", number_max_cnt)
            y_values = rout_temp
            nm2 = trapz(y_values .* conj.(y_values) .* rlist.^2, rlist)
            y_values ./= sqrt.(nm2)
        end
        if sve
            writedlm(fnm*"_$(alph)_n_$(n)_l_$(l)_m_$(m)_.dat", hcat(float(real(rlist)), float(real(y_values .* conj.(y_values)))))
        else
            if return_erg
                return rlist, y_values, erg
            else
                return rlist, y_values
            end
        end
    end
        
end

function radial_inf(erg, mu, M, a, l, m; rpts=1000, rmax_val=1e4, debug=false, iter=50, xtol=1e-120, ftol=1e-120, sve_for_test=false, fnm="test_store/test_radial", is_infin=true)
    # r, erg, unitless
    
    
    rp = 1.0 .+ sqrt.(1.0 .- a.^2)
    r = LinRange(rp .* 1.01, rmax_val .* rp, rpts)
    h = r[2] .- r[1]
    
    alph = mu .* GNew .* M
    k = sqrt.(Complex.(erg.^2 .- alph.^2))
    b = sqrt.(1 - a.^2)
    gam2 = - a * sqrt.(erg.^2 .- alph.^2)
    LLM = l * (l + 1)
    LLM += (-1 + 2 * l * (l + 1) - 2 * m.^2) * gam2 ./ (-3 + 4 * l * (l + 1))
    LLM += ((l - m - 1 * (l - m) * (l + m) * (l + m - 1)) ./ ((-3 + 2 * l) * (2 * l - 1).^2) - (l + 1 - m) * (2 * l - m) * (l + m + 1) * (2 + l + m) ./ ((3 + 2 * l).^2 * (5 + 2 * l))) * gam2.^2 ./ (2 * (1 + 2 * l))
    LLM += (4 * ((-1 + 4 * m^2) * (l * (1 + l) * (121 + l * (1 + l) * (213 + 8 * l * (1 + l) * (-37 + 10 * l * (1 + l)))) - 2 * l * (1 + l) * (-137 + 56 * l * (1 + l) * (3 + 2 * l * (1 + l))) * m^2 + (705 + 8 * l * (1 + l) * (125 + 18 * l * (1 + l))) * m^4 - 15 * (1 + 46 * m^2))) * gam2^3) / ((-5 + 2 * l) * (-3 + 2 * l) * (5 + 2 * l) * (7 + 2 * l) * (-3 + 4 * l * (1 + l))^5)
    
    
    delt = r.^2 .- 2 .* r .+ a.^2
    dr_delt = 2 .* r .- 2
    # rhs = Float64.((erg.^2 .* (r.^2 .+ a.^2).^2 .- 4 .* a .* m .* erg .* r .+ m.^2 .* a.^2) ./ delt .- (erg.^2 .* a.^2 .+ alph.^2 .* r.^2 .+ LLM))
    # rhs = real((erg.^2 .* (r.^2 .+ a.^2).^2 .- 4 .* a .* m .* erg .* r .+ m.^2 .* a.^2) ./ delt .- (erg.^2 .* a.^2 .+ alph.^2 .* r.^2 .+ LLM))
    rhs = ((erg.^2 .* (r.^2 .+ a.^2).^2 .- 4 .* a .* m .* erg .* r .+ m.^2 .* a.^2) ./ delt .- (erg.^2 .* a.^2 .+ alph.^2 .* r.^2 .+ LLM))
    
    # Create the second derivative matrix (Laplacian matrix), and first derivative
    D2 = zeros(rpts, rpts)

    # Fill the diagonal
    for i in 2:rpts-1
        D2[i, i-1] += 1.0 .* delt[i] ./ h^2
        D2[i, i] += -2.0 .* delt[i] ./ h^2
        D2[i, i+1] += 1.0 .* delt[i] ./ h^2
        
        D2[i, i+1] += 1.0 .* dr_delt[i] ./ (2 * h)
        D2[i, i-1] += -1.0 .* dr_delt[i] ./ (2 * h)
    end

    # BNDRY TERMS
    D2[1, 1] += -3.0 .* dr_delt[1] ./ (2 * h)
    D2[1, 2] += 4.0 .* dr_delt[1] ./ (2 * h)
    D2[1, 3] += -1.0 .* dr_delt[1] ./ (2 * h)
    D2[rpts, rpts] += 3.0 .* dr_delt[rpts] ./ (2 * h)
    D2[rpts, rpts - 1] += -4.0 .* dr_delt[rpts] ./ (2 * h)
    D2[rpts, rpts - 2] += 1.0 .* dr_delt[rpts] ./ (2 * h)
    
    D2[1, 1] += 2.0 .* delt[1] ./ h.^3
    D2[1, 2] += -5.0 .* delt[1] ./ h.^3
    D2[1, 3] += 4.0 .* delt[1] ./ h.^3
    D2[1, 4] += -1.0 .* delt[1] ./ h.^3
    D2[rpts, rpts] += 2.0 .* delt[rpts] ./ h.^3
    D2[rpts, rpts - 1] += -5.0 .* delt[rpts] ./ h.^3
    D2[rpts, rpts - 2] += 4.0 .* delt[rpts] ./ h.^3
    D2[rpts, rpts - 3] += -1.0 .* delt[rpts] ./ h.^3
    

    
    function wrapper!(F, x)
        # temp = (D2 * x) .+ rhs .* x
        temp = x.^(-1) .* (D2 * x) .+ rhs
        # F .= real(temp) .+ imag(temp)
        F .= abs.(temp)
        
        
        # real_r = x[1:Int(length(x)/2)]
        # imag_r = x[Int(length(x)/2 + 1):end]
        # r_tot = real_r .+ im .* imag_r
        # temp = (D2 * r_tot) .+ rhs .* r_tot
        # F[1:Int(length(x)/2)] .= real(temp) .* 1e60
        # F[Int(length(x)/2 + 1):end] .= imag(temp) .* 1e60
    end
    
    k = sqrt.(erg.^2 .- alph.^2)
    
    rGuess = radial_inf_NR(k, l, mu, M, r)
    
  
    r_in = real(rGuess)
    
    sol = nlsolve(wrapper!, r_in, show_trace=debug, autodiff = :forward, xtol=xtol, ftol=ftol, iterations=iter)
    full_out = sol.zero
    
   
    if is_infin
        ## need to normalize
        done_search=false
        found1=false
        found2=false
        cnt1 = 0
        cnt2 = 0
        
        indx = rpts
        while !done_search
            if sign(real(full_out[indx])) != sign(real(full_out[indx-1]))
                cnt1 += 1
            end
            if sign(real(r_in[indx])) != sign(real(r_in[indx-1]))
                cnt2 += 1
            end
            
            if cnt1 > 4
                found1 = true
            end
            if cnt2 > 4
                found2 = true
            end
                
            indx -= 1
            if found1&&found2
                done_search = true
            end
        end
        nm1 = maximum(abs.(full_out[indx:end]))
        nm2 = maximum(abs.(r_in[indx:end]))
        
        full_out .*= nm2 ./ nm1
        
        if sve_for_test
            writedlm(fnm*"_real.dat", hcat(r, real(full_out)))
            writedlm(fnm*"_imag.dat", hcat(r, imag(full_out)))
        end
        
        return r, full_out
    else
        
    end
end

function spheroidals(l, m, a, erg)
    # pass erg in normalized units
    if abs.(m) > l
        return 0.0
    end

    # Protection: spin-weighted spheroidal harmonics break down for |m| >= 14
    # In that regime, use spherical harmonics as approximation (valid for small rotation a)
    if abs(m) >= 14
        # Return a wrapper function that evaluates spherical harmonics instead
        return (theta, phi) -> sphericalY(l, m, theta, phi)
    end

    Zlm = spin_weighted_spheroidal_harmonic(0, l, m, a .* erg)
    return Zlm
end


function find_im_part(mu, M, a, n, l, m; debug=false, Ntot_force=5000, iter=50000, xtol=1e-50, ftol=1e-70, return_both=false, for_s_rates=true, QNM=false, QNM_E=1.0, erg_Guess=nothing, max_n_qnm=10, prec=300)
    # output divided by alpha! \tau_nlm = 1 / (2 * output * alpha)!
    
    setprecision(BigFloat, prec)
    
    OmegaH = BigFloat(a ./ (2 .* (GNew .* M) .* (1 .+ sqrt.(1 .- a.^2))))
    alph = BigFloat(mu * GNew * M)
    a = BigFloat(a)
    M = BigFloat(M)
    
    if (ergL(n, l, m, mu, M, a) < m .* OmegaH)||(for_s_rates==true)||QNM
        
        
        
        if m == 1
            Ntot = 3000
        elseif (m==2)
            Ntot = 4250
        elseif (m==3)
            Ntot = 10000
        else
            Ntot = 4000
        end
        
        if a > 0.95
            Ntot *= 2.5
            Ntot = Int(Ntot)
        end
        
        if Ntot_force > Ntot
            Ntot = Ntot_force
        end

        
        
        b = sqrt.(1 - a.^2)
        
        
        if !QNM
            if isnothing(erg_Guess)
                SR211_g = sr_rates(n, l, m, alph ./ (GNew * M), M, a)
                
                if SR211_g == 0
                    SR211_g = -1e-10 .* mu
                end
                w0 = (ergL(n, l, m, alph ./ (GNew * M), M, a; full=false) .+ im * SR211_g) .* GNew * M
                # print("test \t ", w0, "\n")
            else
                w0 = erg_Guess
            end
        else
            ### only valid for l= 0,m=0, otherwise better guesses are needed
            w0 = 0.117 .+ 0.004 .* (alph ./ 0.1) .- 0.084 * im .- 0.01 .* (alph ./ 0.1)  # rough guess for n = 1, l=0, m=0
        end
        
        
        
        function wrapper!(F, x)
            # wR = 10 .^ x[1]  # real
            wR = x[1]
            wI = x[2]  # imag
            erg = wR + im * wI
    
            q = - sqrt.(alph.^2 .- erg.^2)
            
            if !QNM
                if real(q) > 0
                    q *= -1
                end
            else
                if real(q) < 0
                    q *= -1
                end
            end
           
            
            gam = im * a * sqrt.(erg.^2 .- alph.^2)
            LLM = l * (l + 1)
            LLM += (-1 + 2 * l * (l + 1) - 2 * m.^2) * gam.^2 ./ (-3 + 4 * l * (l + 1))
            LLM += ((l - m - 1 * (l - m) * (l + m) * (l + m - 1)) ./ ((-3 + 2 * l) * (2 * l - 1).^2) - (l + 1 - m) * (2 * l - m) * (l + m + 1) * (2 + l + m) ./ ((3 + 2 * l).^2 * (5 + 2 * l))) * gam.^4 ./ (2 * (1 + 2 * l))
            LLM += (4 * ((-1 + 4 * m^2) * (l * (1 + l) * (121 + l * (1 + l) * (213 + 8 * l * (1 + l) * (-37 + 10 * l * (1 + l)))) - 2 * l * (1 + l) * (-137 + 56 * l * (1 + l) * (3 + 2 * l * (1 + l))) * m^2 + (705 + 8 * l * (1 + l) * (125 + 18 * l * (1 + l))) * m^4 - 15 * (1 + 46 * m^2))) * gam^6) / ((-5 + 2 * l) * (-3 + 2 * l) * (5 + 2 * l) * (7 + 2 * l) * (-3 + 4 * l * (1 + l))^5)

            
            c0 = 1.0 .- 2.0 * im * erg - 2 * im ./ b .* (erg .- a .* m ./ 2.0)
            c1 = -4.0 .+ 4 * im * (erg - im * q * (1.0 + b)) + 4 * im / b * (erg - a * m / 2.0) .- 2.0 * (erg.^2 .+ q.^2) ./ q
            c2 = 3.0 - 2 * im * erg - 2.0 * (q.^2 - erg.^2) ./ q - 2.0 * im / b * (erg - a * m ./ 2)
            c3 = 2.0 * im * (erg - im * q).^3 ./ q .+ 2 * (erg .- im * q).^2 .* b + q.^2 * a.^2 .+ 2 * im * q * a * m - LLM - 1 - (erg - im * q).^2 ./ q .+ 2 * q * b + 2 * im / b * ( (erg - im * q).^2 ./ q + 1.0) * (erg .- a * m / 2)
            c4 = (erg .- im * q).^4 ./ q.^2 .+ 2 * im * erg * (erg .- im * q).^2 ./ q .- 2 * im ./ b .* (erg .- im * q).^2 ./ q .* (erg .- a * m ./ 2)
           
            function alphaN(nn)
                return nn.^2 .+ (c0 .+ 1) * nn + c0
            end
            function betaN(nn)
                return -2 * nn.^2 .+ (c1 .+ 2) * nn + c3
            end
            function gammaN(nn)
                return nn.^2 .+ (c2 .- 3) * nn + c4
            end
            
            temp = 1
            for i in Ntot:-1:1
                temp = gammaN(i) ./ (betaN(i) .- alphaN(i) .* temp)
            end
            recur = betaN(0) ./ alphaN(0) .- temp
            
            F[1] = real(recur)
            F[2] = imag(recur)
            
            # print(x, "\t", F, "\n")
        end
        
        if real(w0) < 0
            w0 = alph + imag(w0)
        end
        
        # sol = nlsolve(wrapper!, [BigFloat(log10.(real(w0))), BigFloat(imag(w0))], autodiff = :forward, xtol=xtol, ftol=ftol, iterations=iter)
        if !QNM
            sol = nlsolve(wrapper!, [BigFloat((real(w0))), BigFloat(imag(w0))], autodiff = :forward, xtol=xtol, ftol=ftol, iterations=iter)
        else
            sol = nlsolve(wrapper!, [BigFloat((real(w0))), BigFloat(imag(w0))], autodiff = :forward, xtol=xtol, ftol=ftol, iterations=iter)
            outSend = [n sol.zero[1] sol.zero[2]] # n, real, im
            
            nmax = n + max_n_qnm
            nn_in = n + 1
            while nn_in < nmax
                guess = outSend[end, 2] + (outSend[end, 3] .- 0.1) * im
                if real(guess) < 0.05
                    guess += 0.05
                end
                found_it = false
                trialF = 0
                while !found_it
                    sol = nlsolve(wrapper!, [BigFloat((real(guess))), BigFloat(imag(guess))], autodiff = :forward, xtol=(xtol ./ 1e10) , ftol=ftol, iterations=iter)
                    # print("guess \t", guess, "\t", trialF, "\n")
                    # print("trial \t ", sol.zero, "\n\n")
                    if (!isapprox(sol.zero[1], outSend[end, 2]; atol = 0.0001)&&!isapprox(sol.zero[2], outSend[end, 3]; atol = 0.0001))&&(sol.zero[1] > 0)&&(sol.zero[2] < outSend[end, 3])&&(sol.zero[1] > 1e-3)
                        outSend = cat(outSend, [nn_in sol.zero[1] sol.zero[2]], dims=1)
                        found_it = true
                        if debug
                            print(sol, "\n")
                        end
                    else
                        guess = real(guess) + im .* imag(guess) - 0.005 * im
                    end
                    
                    if (sol.zero[1] == 0)&&(sol.zero[2] == 0)
                       found_it = true
                       nn_in = nmax
                    end
                    
                    trialF +=1
                    if trialF > 200
                        found_it = true
                        nn_in = nmax
                    end
                end
                
                nn_in += 1
            end
            
            return outSend
        end
        
        if debug
            print(sol, "\n\n")

          
            Ff = zeros(2)
            wrapper!(Ff, sol.zero)
            print("sanity check \t", Ff, "\n")
            if (abs.(Ff[1]) .> 1e-30)||(abs.(Ff[2]) .> 1e-30)
                if return_both
                    return 0.0, 0.0
                else
                    return 0.0
                end
            end
        end
        
        if ((!sol.x_converged) && (!sol.f_converged))||(sol.zero[1] .< 0.0)
            if return_both
                return 0.0, 0.0
            else
                return 0.0
            end
        end
        
        

        if !return_both
            return sol.zero[2]   # im(G * M * omega)
        else
            return sol.zero[1], sol.zero[2]
        end
    else
        if !return_both
            return 0.0
        else
            return 0.0, 0.0
        end
    end
end


function find_im_zero(mu, M, n, l, m; debug=false, Ntot=3000, iter=1000, xtol=1e-6, ftol=1e-20)
    
    
    alph = BigFloat(mu * GNew * M)
    QNM = false
        
    function wrapper!(F, x)
        # wR = 10 .^ x[1]  # real
        wR = x[1]
        a =  x[2]
        
        
        wI = 0.0
        erg = wR + im * wI
        
        b = sqrt.(Complex(1 - a.^2))

        q = - sqrt.(Complex(alph.^2 .- erg.^2))
        
        if !QNM
            if real(q) > 0
                q *= -1
            end
        else
            if real(q) < 0
                q *= -1
            end
        end
       
        
        gam = im * a * sqrt.(Complex(erg.^2 .- alph.^2))
        LLM = l * (l + 1)
        LLM += (-1 + 2 * l * (l + 1) - 2 * m.^2) * gam.^2 ./ (-3 + 4 * l * (l + 1))
        LLM += ((l - m - 1 * (l - m) * (l + m) * (l + m - 1)) ./ ((-3 + 2 * l) * (2 * l - 1).^2) - (l + 1 - m) * (2 * l - m) * (l + m + 1) * (2 + l + m) ./ ((3 + 2 * l).^2 * (5 + 2 * l))) * gam.^4 ./ (2 * (1 + 2 * l))
        LLM += (4 * ((-1 + 4 * m^2) * (l * (1 + l) * (121 + l * (1 + l) * (213 + 8 * l * (1 + l) * (-37 + 10 * l * (1 + l)))) - 2 * l * (1 + l) * (-137 + 56 * l * (1 + l) * (3 + 2 * l * (1 + l))) * m^2 + (705 + 8 * l * (1 + l) * (125 + 18 * l * (1 + l))) * m^4 - 15 * (1 + 46 * m^2))) * gam^6) / ((-5 + 2 * l) * (-3 + 2 * l) * (5 + 2 * l) * (7 + 2 * l) * (-3 + 4 * l * (1 + l))^5)

        
        c0 = 1.0 .- 2.0 * im * erg - 2 * im ./ b .* (erg .- a .* m ./ 2.0)
        c1 = -4.0 .+ 4 * im * (erg - im * q * (1.0 + b)) + 4 * im / b * (erg - a * m / 2.0) .- 2.0 * (erg.^2 .+ q.^2) ./ q
        c2 = 3.0 - 2 * im * erg - 2.0 * (q.^2 - erg.^2) ./ q - 2.0 * im / b * (erg - a * m ./ 2)
        c3 = 2.0 * im * (erg - im * q).^3 ./ q .+ 2 * (erg .- im * q).^2 .* b + q.^2 * a.^2 .+ 2 * im * q * a * m - LLM - 1 - (erg - im * q).^2 ./ q .+ 2 * q * b + 2 * im / b * ( (erg - im * q).^2 ./ q + 1.0) * (erg .- a * m / 2)
        c4 = (erg .- im * q).^4 ./ q.^2 .+ 2 * im * erg * (erg .- im * q).^2 ./ q .- 2 * im ./ b .* (erg .- im * q).^2 ./ q .* (erg .- a * m ./ 2)
       
        function alphaN(nn)
            return nn.^2 .+ (c0 .+ 1) * nn + c0
        end
        function betaN(nn)
            return -2 * nn.^2 .+ (c1 .+ 2) * nn + c3
        end
        function gammaN(nn)
            return nn.^2 .+ (c2 .- 3) * nn + c4
        end
        
        temp = 1
        for i in Ntot:-1:1
            temp = gammaN(i) ./ (betaN(i) .- alphaN(i) .* temp)
        end
        recur = betaN(0) ./ alphaN(0) .- temp
        
        F[1] = real(recur)
        F[2] = imag(recur)
        
        # print(x, "\t", F, "\n")
    end
        
        
    w0G = alph .* (1 .- alph.^2 ./ n.^2)
    sol = nlsolve(wrapper!, [w0G, 0.1], autodiff = :forward, xtol=xtol, ftol=ftol, iterations=iter, show_trace=debug)
    println(sol)
    return sol.zero[2]
       
end



function compute_gridded(mu, M, a, n, l, m; Ntot=2000, iter=50, xtol=1e-7, npts=30, amin=0.0, compute_neg=false)
    a_max = a * 1.1 # safety, just in case upward fluctuation
    if a_max > maxSpin
        a_max = maxSpin
    end
    
    output = zeros(npts)
    if !compute_neg
        if amin == a_max
            amin *= 0.99
        end
        if amin < a_max
            alist = LinRange(amin, a_max, npts);
        
            alph = GNew * M * mu

            for i in 1:length(alist)
                output[i] = find_im_part(mu, M, alist[i], n, l, m, Ntot_force=Ntot, iter=iter, xtol=xtol, for_s_rates=true) ./ (GNew * M)
            end
        else
            alist = LinRange(a_max, amin, npts);
        end
    else
        if amin > maxSpin
            amin = maxSpin
        end
        alist = LinRange(minimum([0.1 amin]), maximum([0.1 amin]), npts);
        
        alph = GNew * M * mu

        for i in 1:length(alist)
            output[i] = find_im_part(mu, M, alist[i], n, l, m, Ntot_force=Ntot, iter=iter, xtol=xtol, for_s_rates=true) ./ (GNew * M)
        end
    end
    return alist, output
end

function generalized_laguerre(n, α, x)
    return (-1).^ n ./ factorial(big(n)) .* HypergeometricFunctions.U.(-n, α+1 ,x)
    # return binomial(n + α, n) .* HypergeometricFunctions.M.(-n, α+1 , x)
end

function test_projection_scatter(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=20000, Npts_Bnd=4000, rmaxT=100, debug=false, Ntot_safe=5000, sve_for_test=false)
    
    rp = 1 + sqrt.(1 - a.^2)
    alph = mu * GNew * M
    
    # approx energy estimates
    # erg_1G = ergL(n1, l1, m1, mu, M, a)
    # erg_2G = ergL(n2, l2, m2, mu, M, a)
    # erg_3G = ergL(n3, l3, m3, mu, M, a)
    # erg_pxy = (erg_1G + erg_2G - erg_3G) * GNew * M
    
    wR1, wI1 = find_im_part(mu, M, a, n1, l1, m1; return_both=true, for_s_rates=true)
    wR2, wI2 = find_im_part(mu, M, a, n2, l2, m2; return_both=true, for_s_rates=true)
    wR3, wI3 = find_im_part(mu, M, a, n3, l3, m3; return_both=true, for_s_rates=true)
    erg_pxy = (wR1 + wR2 - wR3)
    
    qq = abs.(real(- sqrt.(alph.^2 .- erg_pxy.^2)))
    rmax = rmaxT ./ qq
    
    # rlist = 10 .^(range(log10.(rp), log10.(rmax), rpts))
    rlist = range(rp, rmax, rpts)
    
    theta_rand = 2 * acos.(sqrt.(rand(rpts)))
    # r_rand = rand(rpts).^(1/3) * (rmax^3 ./ 3 - 1.0 / 3)^(1/3)
    sampls = rand(rpts)
    r_rand = rp .+ sampls .* (rmax .- rp)
    exF = (rmax .- rp)
    
    alph = mu .* GNew .* M
    
    
    
    erg_out = ergL(1, 0, 0, mu, M, a)
    erg_diff = abs.(-(erg_pxy.^2 .- (erg_out .* GNew .* M).^2))
    
    idx_run = 1
    for i in 2:30
        erg_out = ergL(i, 0, 0, mu, M, a)
        erg_diff_hold = abs.(-(erg_pxy.^2 .- (erg_out .* GNew .* M).^2))
        if erg_diff_hold < erg_diff
            erg_diff = erg_diff_hold
            idx_run = i
        end
    end
    
    
    for i in idx_run:idx_run
        rl, r1, erg = solve_radial(mu, M, a, i, 0, 0; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
        
        itp_r = LinearInterpolation(log10.(Float64.(rl)), real(r1), extrapolation_bc=Line())
        itp_im = LinearInterpolation(log10.(Float64.(rl)), imag(r1), extrapolation_bc=Line())
        function rad_comp(r)
            return itp_r(log10.(r)) .+ im * itp_im(log10.(r))
        end
        

        # erg = ergL(i, 0, 0, mu, M, a) .* GNew .* M
        erg = erg_pxy
        
        r = LinRange(rp .* 1.01, rmax, rpts)
        h = r[2] .- r[1]
        
        alph = mu .* GNew .* M
        k = sqrt.(Complex.(erg.^2 .- alph.^2))
                
        delt = r.^2 .- 2 .* r .+ a.^2
        dr_delt = 2 .* r .- 2
        m = 0
        # rhs = real((erg.^2 .* (r.^2 .+ a.^2).^2 .- 4 .* a .* m .* erg .* r .+ m.^2 .* a.^2) ./ delt .- (erg.^2 .* a.^2 .+ alph.^2 .* r.^2))
        rhs = ((erg.^2 .* (r.^2 .+ a.^2).^2 .- 4 .* a .* m .* erg .* r .+ m.^2 .* a.^2) ./ delt .- (erg.^2 .* a.^2 .+ alph.^2 .* r.^2))
        
        # Create the second derivative matrix (Laplacian matrix), and first derivative
        D2 = zeros(rpts, rpts)

        # Fill the diagonal
        for i in 2:rpts-1
            D2[i, i-1] += 1.0 .* delt[i] ./ h^2
            D2[i, i] += -2.0 .* delt[i] ./ h^2
            D2[i, i+1] += 1.0 .* delt[i] ./ h^2
            
            D2[i, i+1] += 1.0 .* dr_delt[i] ./ (2 * h)
            D2[i, i-1] += -1.0 .* dr_delt[i] ./ (2 * h)
        end

        # BNDRY TERMS
        D2[1, 1] += -1.0 .* dr_delt[1] ./ (2 * h)
        D2[1, 2] += 1.0 .* dr_delt[1] ./ (2 * h)
        # D2[1, 3] += -1.0 .* dr_delt[1] ./ (2 * h)
        D2[rpts, rpts] += 1.0 .* dr_delt[rpts] ./ (2 * h)
        D2[rpts, rpts - 1] += -1.0 .* dr_delt[rpts] ./ (2 * h)
        # D2[rpts, rpts - 2] += 1.0 .* dr_delt[rpts] ./ (2 * h)
        
        D2[1, 1] += 1.0 .* delt[1] ./ h.^2
        D2[1, 2] += -2.0 .* delt[1] ./ h.^2
        D2[1, 3] += 1.0 .* delt[1] ./ h.^2
        # D2[1, 4] += -1.0 .* delt[1] ./ h.^2
        D2[rpts, rpts] += 1.0 .* delt[rpts] ./ h.^2
        D2[rpts, rpts - 1] += -2.0 .* delt[rpts] ./ h.^2
        D2[rpts, rpts - 2] += 1.0 .* delt[rpts] ./ h.^2
        # D2[rpts, rpts - 3] += -1.0 .* delt[rpts] ./ h.^2
        
        r_TEST = rad_comp(r)
        temp = (D2 * r_TEST) .+ rhs .* r_TEST
      
        out = Float64(real(-trapz(temp .* conj(r_TEST) , r)))
 
        ################
        
        
        # writedlm("test_store/check_real_1.dat", hcat(real(rlist), float(real(itp_dr(log10.(rlist))))))
        # writedlm("test_store/check_im_1.dat", hcat(real(rlist), float(imag(itp_dr(log10.(rlist))))))
        # writedlm("test_store/check_real_2.dat", hcat(real(rlist), float(real(itp_dr2(log10.(rlist))))))
        # writedlm("test_store/check_im_2.dat", hcat(real(rlist), float(imag(itp_dr2(log10.(rlist))))))
        

        # print(integrd, "\n")
        if sve_for_test
            writedlm("test_store/check_real_1.dat", hcat(real(r_rand), float(real(d_t2phi .+ dr_term .+ theta_term .+ mu_term .* Z_outC .* conj(rad_comp(r_rand))))))
        end

        
        erg_out = ergL(i, 0, 0, mu, M, a)
        erg_diff_pxy = -(erg_pxy.^2 .- (erg_out .* GNew .* M).^2)
        
        
        
#        if debug
#            if i == 1
#                print(i, "\t", -11 * alph^4 ./ 18, "\t", erg_diff_pxy, "\t", out, "\n")
#            elseif i == 2
#                print(i, "\t", 5 * alph^4 ./ 36, "\t", erg_diff_pxy, "\t", out, "\n")
#            end
#        end


        return erg_diff_pxy, out
        
    end
end



function gf_radial(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=1000, Npts_Bnd=1000, debug=false, Ntot_safe=10000,  iter=10, xtol=1e-10, ftol=1e-10, tag="_", Nang=500000, eps_fac=1e-10, NON_REL=false, h_mve=1, to_inf=false, rmaxT=100, prec=200, cvg_acc=1e-4, NptsCh=60, iterC=40, Lcheb=4, run_leaver=false, der_acc=1e-20, use_heunc=false)

    
    rp = BigFloat(1.0 .+ sqrt.(1.0 .- a.^2))
    rmm = BigFloat(1.0 .- sqrt.(1.0 .- a.^2))
    alph = mu * GNew * M
    
    
    
    
    maxN = maximum([n1 n2 n3])
    minN = maximum([n1 n2 n3])
    
    erg_pxy = 0.0
    erg_1G = 0.0
    erg_2G = 0.0
    erg_3G = 0.0
    ### quick check if to inf to not!
    try
        erg_1G, erg1I = find_im_part(mu, M, a, n1, l1, m1; Ntot_force=Ntot_safe, for_s_rates=true, return_both=true)
        erg_2G, erg2I = find_im_part(mu, M, a, n2, l2, m2; Ntot_force=Ntot_safe, for_s_rates=true, return_both=true)
        erg_3G, erg3I = find_im_part(mu, M, a, n3, l3, m3; Ntot_force=Ntot_safe, for_s_rates=true, return_both=true)
        
        if erg_1G == 0 || erg_2G == 0 || erg_3G == 0
            erg_1G = ergL(n1, l1, m1, mu, M, a; full=true) .* GNew .* M
            erg_2G = ergL(n1, l1, m1, mu, M, a; full=true) .* GNew .* M
            erg_3G = ergL(n1, l1, m1, mu, M, a; full=true) .* GNew .* M
        end
        
        erg_pxy = (erg_1G + erg_2G - erg_3G)
    catch
        erg_1G = ergL(n1, l1, m1, mu, M, a; full=true)
        erg_2G = ergL(n2, l2, m2, mu, M, a; full=true)
        erg_3G = ergL(n3, l3, m3, mu, M, a; full=true)
        erg_pxy = (erg_1G + erg_2G - erg_3G) .* GNew .* M
    end
    
    OmegaH = a ./ (2 .* (GNew .* M) .* (1 .+ sqrt.(1 .- a.^2)))
    
    if erg_pxy .> alph
        to_inf = true
    else
        to_inf = false
    end
    
    if debug
        println("erg fin proxy \t", erg_pxy)
    end
    
    #### fix l,m
    if to_inf
        m = (m1 + m2 - m3)
        l = l1 + l2 - l3
        
        while !((l >= m)&&(l>0))
            l += 2
        end
        
        
        kk_pxy = real(sqrt.(erg_pxy.^2 .- alph.^2))
        rmax = 1/kk_pxy .* 20.0 # Need prefactor here i think...
        
        val, iszero = integral4(l1, m1, l2, m2, l3, -m3, l, -m)
        tcnt = 0
        while iszero
            val, iszero = integral4(l1, m1, l2, m2, l3, -m3, l, -m)
            if iszero
                l += 1
                tcnt += 1
                if tcnt > 10
                    iszero = false
                    return 0.0
                end
            end
        end

    else
        l = 0
        m = 0
        # rmax = Float64.(100 ./ alph.^2 .* (minN ./ 2.0) )
        if maxN < 10
            rmax = 2.0 .^(2.0 .* maxN .- 2 .* (1 .+ maxN)) .* gamma(2 .+ 2 .* maxN) ./ alph.^2 ./ factorial(2 .* maxN - 1) .* 10.0
        else
            rmax = 2.0 .^(2.0 .* maxN .- 2 .* (1 .+ maxN)) .* 4 .* maxN .^2 ./ alph.^2 .* 10.0
        end
    end
    
    
    rlist = 10 .^range(log10(rp .* (1.0 .+ eps_fac)), log10.(rmax), rpts)
    
    # simplify_radial = check_slv_rad(alph, m1)
    simplify_radial = false
    if NON_REL||simplify_radial
        rf_1 = radial_bound_NR(n1, l1, m1, mu, M, rlist)
        erg_1 = erg_1G
    else
        if run_leaver
            rl, r1, erg_1 = solve_radial(mu, M, a, n1, l1, m1; rpts=Npts_Bnd, return_erg=true, Ntot_safe=Ntot_safe, use_heunc=use_heunc, debug=debug)
        else
            wR, wI, rl, r1 = eigensys_Cheby(M, a, mu, n1, l1, m1, return_wf=true, Npoints=NptsCh, Iter=iterC, L=Lcheb, cvg_acc=cvg_acc, Npts_r=Npts_Bnd, return_nu=false, prec=prec, sfty_run=false, der_acc=der_acc, debug=debug)
            erg_1 = wR .+ im .* wI
        end
        itp = LinearInterpolation(log10.(rl), r1, extrapolation_bc=Line())
        rf_1 = itp(log10.(rlist))
        real_diff_test = abs.((erg_1 .- erg_1G) ./ alph) # saftey net for random fail....
        if (imag(erg_1) < 0)||(real_diff_test .> 0.5)
            if debug
                println("Issue with one of the energy eigenstates....")
                println("ERG 1\t", erg_1)
            end
            if (erg_1G .< 0.8 .* m1 .* OmegaH) # if fail is deep in non-relativistic regime, catch
                rf_1 = radial_bound_NR(n1, l1, m1, mu, M, rlist)
                erg_1 = erg_1G
            else
                return 0.0
            end
        end
    end
    
    
    
    # simplify_radial = check_slv_rad(alph, m2)
    simplify_radial = false
    if NON_REL||simplify_radial
        rf_2 = radial_bound_NR(n2, l2, m2, mu, M, rlist)
        erg_2 = erg_2G
    else
        if (n2 == n1)&&(l2 == l1)&&(m2 == m1)
            rf_2 = rf_1
            erg_2 = erg_1
        else
            if run_leaver
                rl, r2, erg_2 = solve_radial(mu, M, a, n2, l2, m2; rpts=Npts_Bnd,  return_erg=true, Ntot_safe=Ntot_safe, use_heunc=use_heunc, debug=debug)
            else
                wR, wI, rl, r2 = eigensys_Cheby(M, a, mu, n2, l2, m2, return_wf=true, Npoints=NptsCh, Iter=iterC, cvg_acc=cvg_acc, L=Lcheb,  Npts_r=Npts_Bnd, return_nu=false, prec=prec, sfty_run=false, der_acc=der_acc, debug=debug)
                erg_2 = wR .+ im .* wI
            end
            
            itp = LinearInterpolation(log10.(rl), r2, extrapolation_bc=Line())
            rf_2 = itp(log10.(rlist))
            real_diff_test = abs.((erg_2 .- erg_2G ) ./ alph) # saftey net for random fail....
            if (imag(erg_2) < 0)||(real_diff_test .> 0.5)
                if debug
                    println("Issue with one of the energy eigenstates....")
                    println("ERG 2\t", erg_2)
                end
                if (erg_2G .< 0.8 .* m2 .* OmegaH) # if fail is deep in non-relativistic regime, catch
                    rf_2 = radial_bound_NR(n2, l2, m2, mu, M, rlist)
                    erg_2 = erg_2G
                else
                    return 0.0
                end
#                return 0.0
            end
        end
    end
    
    # simplify_radial = check_slv_rad(alph, m3)
    simplify_radial = false
    if NON_REL||simplify_radial
        rf_3 = radial_bound_NR(n3, l3, m3, mu, M, rlist)
        erg_3 = erg_3G
    else
        if run_leaver
            rl, r3, erg_3 = solve_radial(mu, M, a, n3, l3, m3; rpts=Npts_Bnd, return_erg=true, Ntot_safe=Ntot_safe, use_heunc=use_heunc, debug=debug)
        else
            wR, wI, rl, r3 = eigensys_Cheby(M, a, mu, n3, l3, m3, return_wf=true, Npoints=NptsCh, Iter=iterC, cvg_acc=cvg_acc, L=Lcheb, Npts_r=Npts_Bnd, return_nu=false, prec=prec, sfty_run=false, der_acc=der_acc, debug=debug)
            erg_3 = wR .+ im .* wI
        end
        
        itp = LinearInterpolation(log10.(rl), r3, extrapolation_bc=Line())
        rf_3 = itp(log10.(rlist))
        real_diff_test = abs.((erg_3 .- erg_3G) ./ alph) # saftey net for random fail....
        if (imag(erg_3) < 0)||(real_diff_test .> 0.5)
            if debug
                println("Issue with one of the energy eigenstates....")
                println("ERG 3\t", erg_3)
            end
            if (erg_3G .< 0.8 .* m3 .* OmegaH) # if fail is deep in non-relativistic regime, catch
                rf_3 = radial_bound_NR(n3, l3, m3, mu, M, rlist)
                erg_3 = erg_3G
            else
                return 0.0
            end
#            return 0.0
        end
    end
    
    
    # writedlm("test_store/test1.dat", hcat(real(rlist), real(rf_1 .* conj.(rf_1))))
    # writedlm("test_store/test2.dat", hcat(real(rlist), real(rf_2 .* conj.(rf_2))))
    # writedlm("test_store/test3.dat", hcat(real(rlist), real(rf_3 .* conj.(rf_3))))
    erg = (erg_1 + erg_2 - erg_3) + 0 * im # leave the 0 im for NR case
    
    ### recheck if to inf holds...
    if (abs.(erg) .> alph)&&!to_inf
        println("Flipping from bnd to inf")
        to_inf = true
        m = (m1 + m2 - m3)
        l = l1 + l2 - l3
        while l < m
            l += 2
        end
        val, iszero = integral4(l1, m1, l2, m2, l3, -m3, l, -m)
        if iszero
            l += 2
        end
    elseif (abs.(erg) .< alph)&&to_inf
        println("Flipping from inf to bnd")
        to_inf = false
        l = 0
        m = 0
    end
    
    Z1 = spheroidals(l1, m1, a, erg_1)
    Z2 = spheroidals(l2, m2, a, erg_2)
    Z3 = spheroidals(l3, m3, a, erg_3)
    Z4 = spheroidals(l, m, a, erg)

    function compute_angular_integral_adaptive(Z1, Z2, Z3, Z4, use_spherical_fallback=false, timeout_seconds=30.0)
        """
        Adaptive Monte Carlo integration with timeout and fallback to spherical harmonics.
        Returns (CG, CG_2, used_fallback) where CG and CG_2 are the integrated values.
        """

        # Simple cache for spheroidal evaluations (key = (theta, phi) rounded to grid, value = product)
        cache = Dict{Tuple{Float64, Float64}, Tuple{Float64, Float64}}()
        CACHE_GRID = 1000  # Round to ~3 decimal places for caching

        function cache_key(theta::Float64, phi::Float64)
            return (round(theta * CACHE_GRID) / CACHE_GRID, round(phi * CACHE_GRID) / CACHE_GRID)
        end

        # Precompute spheroidal products to avoid redundant evaluations
        function func_ang_pair(theta, phi)
            try
                key = cache_key(theta, phi)
                if haskey(cache, key)
                    return cache[key]
                end

                if use_spherical_fallback
                    # Use spherical harmonics (much faster)
                    S1 = sphericalY(l1, m1, theta, phi)
                    S2 = sphericalY(l2, m2, theta, phi)
                    S3 = sphericalY(l3, m3, theta, phi)
                    S4 = sphericalY(l, m, theta, phi)
                    Z_prod = S1 * S2 * conj(S3) * conj(S4)
                else
                    Z_prod = Z1(theta, phi) * Z2(theta, phi) * conj(Z3(theta, phi)) * conj(Z4(theta, phi))
                end
                real_prod = real(Z_prod)
                result = (real_prod, real_prod * cos(theta)^2)
                cache[key] = result
                return result
            catch
                return 0.0, 0.0
            end
        end

        CG = 0.0
        CG_2 = 0.0
        ang_cvrg = false
        idx = 0
        check_interval = 100  # Check convergence every 100 samples
        CG_old = 0.0
        CG_2_old = 0.0
        max_iterations = 100000  # Safeguard against infinite loops
        tol = 1e-4
        min_samples = 200  # Require minimum samples before checking convergence
        max_samples = 50000  # Stop if we exceed this many samples

        # Track time to implement timeout
        start_time = time()

        while !ang_cvrg && idx < max_iterations
            # Check timeout periodically
            if idx % 1000 == 0 && (time() - start_time) > timeout_seconds
                if debug
                    println("Angular MC integration timeout! Using spherical harmonics fallback.")
                end
                # Fallback to spherical harmonics
                if !use_spherical_fallback
                    return compute_angular_integral_adaptive(Z1, Z2, Z3, Z4, true, timeout_seconds)
                else
                    # Already tried fallback, return best estimate with warning
                    if idx > 0
                        CG *= 4*pi / idx
                        CG_2 *= 4*pi / idx
                    end
                    return CG, CG_2, true
                end
            end

            # Batch evaluate multiple samples per iteration
            for _ in 1:check_interval
                theta = acos(1.0 - 2.0 * rand())
                phi = 2.0 * pi * rand()
                f1, f2 = func_ang_pair(theta, phi)
                CG += f1
                CG_2 += f2
                idx += 1
            end

            if idx >= min_samples
                # Compute running averages
                CG_avg = CG / idx
                CG_2_avg = CG_2 / idx

                # Check convergence with safeguard against division by zero
                if abs(CG_old) > 1e-10
                    test1 = abs(CG_avg - CG_old) / abs(CG_old)
                else
                    test1 = abs(CG_avg - CG_old)
                end

                if abs(CG_2_old) > 1e-10
                    test2 = abs(CG_2_avg - CG_2_old) / abs(CG_2_old)
                else
                    test2 = abs(CG_2_avg - CG_2_old)
                end

                if (test1 < tol) && (test2 < tol)
                    ang_cvrg = true
                elseif idx > max_samples
                    # Force convergence if we've reached max samples
                    ang_cvrg = true
                else
                    CG_old = CG_avg
                    CG_2_old = CG_2_avg
                end
            end
        end

        # Normalize by volume of integration domain (surface of unit sphere = 4π)
        if idx > 0
            CG *= 4*pi / idx
            CG_2 *= 4*pi / idx
        end

        return CG, CG_2, use_spherical_fallback
    end

    println("getting angles")
    CG, CG_2, used_fallback = compute_angular_integral_adaptive(Z1, Z2, Z3, Z4, false, 30.0)

    if used_fallback && debug
        println("Warning: Angular integration used spherical harmonic fallback")
    end

    
    
    if debug
        println("ergs ", erg_1, "  ", erg_2, "  ", erg_3, "  ", erg)
        println("angles \t", CG, "\t", CG_2)
    end
    
    lam_eff = 1.0e0
    if (n1==n2)&&(l1==l2)&&(m1==m2)
        preFac = 1.0 ./ 2.0
    else
        preFac = 1
    end
    unitMatch = -1.0 ./ (2 .* alph).^(3/2) .* lam_eff
    gammaT = (preFac .* (rf_1 .* rf_2 .* conj(rf_3)) .* unitMatch)
        
    itpG = LinearInterpolation(log10.(rlist), Float64.(real.(gammaT)), extrapolation_bc=Line())
    itpGI = LinearInterpolation(log10.(rlist), Float64.(imag.(gammaT)), extrapolation_bc=Line())
    
    gam = im * a * sqrt.(erg.^2 .- alph.^2)
    LLM = l * (l + 1)
    LLM += (-1 + 2 * l * (l + 1) - 2 * m.^2) * gam.^2 ./ (-3 + 4 * l * (l + 1))
    LLM += ((l - m - 1 * (l - m) * (l + m) * (l + m - 1)) ./ ((-3 + 2 * l) * (2 * l - 1).^2) - (l + 1 - m) * (2 * l - m) * (l + m + 1) * (2 + l + m) ./ ((3 + 2 * l).^2 * (5 + 2 * l))) * gam.^4 ./ (2 * (1 + 2 * l))
    LLM += (4 * ((-1 + 4 * m^2) * (l * (1 + l) * (121 + l * (1 + l) * (213 + 8 * l * (1 + l) * (-37 + 10 * l * (1 + l)))) - 2 * l * (1 + l) * (-137 + 56 * l * (1 + l) * (3 + 2 * l * (1 + l))) * m^2 + (705 + 8 * l * (1 + l) * (125 + 18 * l * (1 + l))) * m^4 - 15 * (1 + 46 * m^2))) * gam^6) / ((-5 + 2 * l) * (-3 + 2 * l) * (5 + 2 * l) * (7 + 2 * l) * (-3 + 4 * l * (1 + l))^5)
    
    
    r_list_map = 10 .^LinRange(log10.(rp * (1.0 .+ eps_fac)), log10.(rmax), 100000)
    rout_star = r_list_map .+ 2.0 .* rp ./ (rp .- rmm) .* log.((r_list_map .- rp) ./ 2.0) .- 2.0 .* rmm .* log.((r_list_map .- rmm) ./ 2.0) ./ (rp .- rmm)
    
    
    itp_rrstar = LinearInterpolation(rout_star, r_list_map, extrapolation_bc=Line())
    itp_rrstar_inv = LinearInterpolation(r_list_map, rout_star, extrapolation_bc=Line())
    rr = itp_rrstar_inv(rmax)
    
    h_step = Float64.(rp ./ h_mve)

    # Get solution #1
    outWF = []
    rvals = []
    # set D to 1
    
    append!(outWF,  exp.(-sqrt.(alph.^2 .- erg.^2) .* rr) .* sqrt.(itp_rrstar.(rr).^2 .+ a.^2) ./ itp_rrstar.(rr)) ##
    append!(rvals, rr)
    rr -= h_step
    append!(outWF,  exp.(-sqrt.(alph.^2 .- erg.^2) .* rr) .* sqrt.(itp_rrstar.(rr).^2 .+ a.^2) ./ itp_rrstar.(rr))
    append!(rvals, rr)
    rr -= h_step

    
    stop_running = false
    run_it = true
    idx = 2
    while run_it
        r_input = itp_rrstar(rvals[idx])
    
        delt = (r_input.^2 .- 2 .* r_input .+ a.^2)
        ff = (r_input.^2 .+ a.^2)
            
        net_rescale = h_step.^2
        Vv = delt .* alph.^2 ./ ff .+ delt .* (LLM .+ a.^2 .* (erg.^2 .- alph.^2)) ./ ff.^2 .+ delt .* (3 .* r_input.^2 .- 4 .* r_input .+ a.^2) ./ ff.^3 .- 3 .* delt.^2 .* r_input.^2 ./ ff.^4
        newV = 2 * outWF[idx] .- outWF[idx - 1] .+ net_rescale .* (Vv .- erg.^2) .* outWF[idx]
        
        newV_r = Float64.(real(newV))
        newV_i = Float64.(imag(newV))
        append!(outWF, newV_r + im * newV_i)
            
        append!(rvals, rr)

        rr -= h_step
        idx += 1
        
        if (itp_rrstar(rr) < rp .*  (1.0 .+ eps_fac))
            run_it = false
        end
    end
    rvals = reverse(rvals)
    outWF = reverse(outWF)
    
    # Get solution #2
    rr += h_step
    outWF_fw = []
    append!(outWF_fw, exp.(- im .* erg .* rr) .* sqrt.(itp_rrstar.(rr).^2 .+ a.^2)) ##
    rr += h_step
    append!(outWF_fw, exp.(- im .* erg .* rr) .* sqrt.(itp_rrstar.(rr).^2 .+ a.^2))
    rr += h_step
    
    stop_running = false
    run_it = true
    idx = 2
    while run_it
        r_input = itp_rrstar(rvals[idx])
    
        delt = (r_input.^2 .- 2 .* r_input .+ a.^2)
        ff = (r_input.^2 .+ a.^2)
            
        net_rescale = h_step.^2
        Vv = delt .* alph.^2 ./ ff .+ delt .* (LLM .+ a.^2 .* (erg.^2 .- alph.^2)) ./ ff.^2 .+ delt .* (3 .* r_input.^2 .- 4 .* r_input .+ a.^2) ./ ff.^3 .- 3 .* delt.^2 .* r_input.^2 ./ ff.^4
        newV = 2 * outWF_fw[idx] .- outWF_fw[idx - 1] .+ net_rescale .* (Vv .- erg.^2) .* outWF_fw[idx]
        
        newV_r = Float64.(real(newV))
        newV_i = Float64.(imag(newV))
        append!(outWF_fw, newV_r + im * newV_i)
            
        rr += h_step
        idx += 1
        
        if rr > rvals[end]
            run_it = false
        end
    end
    
    
    outWF .*= 1.0 ./ sqrt.(itp_rrstar.(rvals).^2 .+ a.^2)
    outWF_fw .*= 1.0 ./ sqrt.(itp_rrstar.(rvals).^2 .+ a.^2)

    if length(outWF_fw) > 15
        midP = Int(round(length(rvals) - 10))
    elseif length(outWF_fw) > 4
        midP = Int(round(length(rvals) - 3))
    else
        return 0.0
    end
    wronk  = (outWF_fw[midP] .* (outWF[midP+1] .- outWF[midP-1]) .-  outWF[midP] .* (outWF_fw[midP+1] .- outWF_fw[midP-1]) )./ (itp_rrstar.(rvals[midP+1]) .- itp_rrstar.(rvals[midP-1]))
    wronk *= (itp_rrstar.(rvals[midP]).^2 .- 2 .* itp_rrstar.(rvals[midP]) .+ a.^2) .* (GNew .* M)

    Tmm = (itpG(log10.(itp_rrstar.(rvals))) + im * itpGI(log10.(itp_rrstar.(rvals)))) .* (CG .* itp_rrstar.(rvals).^2 .+ CG_2 .* a.^2)
    
    ####
    if to_inf
        
        out_R = outWF[end] .* trapz(outWF_fw .* Tmm, itp_rrstar.(rvals)) ./ wronk
        maxV = Float64.(real(out_R .* conj.(out_R)))

        lam = (mu ./ (M_pl .* 1e9))^2
        kk = real(sqrt.(erg.^2 .- alph.^2))
        println("1/kk \t", 1 ./ kk)
        rate_out = 2 .* alph .* kk .* (maxV[end] .* itp_rrstar.(itp_rrstar.(rvals[end])).^2) .* lam^2
        out_gamma = rate_out ./ mu^2 .* (GNew * M^2 * M_to_eV)^2
    else
        idx_hold = itp_rrstar.(rvals) .> 1.01 .* rp
        rnew_rp = trapz(outWF[idx_hold] .* Tmm[idx_hold] , itp_rrstar.(rvals[idx_hold])) .* outWF_fw[idx_hold][1] ./ wronk

        maxV = real(rnew_rp.* conj.(rnew_rp))

        lam = (mu ./ (M_pl .* 1e9))^2
        rate_out = 4 .* alph.^2 .* (1 .+ sqrt.(1 - a.^2)) .* Float64(maxV) .* lam^2
        out_gamma = rate_out ./ mu^2 .* (GNew * M^2 * M_to_eV)^2
    end

    if debug
        print("Rate \t", out_gamma, "\n")
    end
    return out_gamma # unitless [gamma / mu]
end


function pre_computed_sr_rates(n, l, m, alph, M; n_high=20, n_low=20, delt_a=0.001, cheby=true)
    state_str = format_state_for_filename(n, l, m)
    if !cheby
        fn = "rate_sve/Imag_zero_$(state_str).dat"
    else
        fn = "rate_sve/Imag_zeroC_$(state_str).dat"
    end
    if isfile(fn)
        zerolist= readdlm(fn)
        
        itp = LinearInterpolation(zerolist[:, 1], zerolist[:, 2], extrapolation_bc=Line())
        a_mid = itp(alph)
        
        if a_mid .< 0
            a_mid = 0.01
        elseif a_mid .> maxSpin
            a_mid = maxSpin
        end
    else
        a_mid = 0.5
    end
    
    run_high = true
    run_low = true
    if (a_mid + delt_a) .> maxSpin
        run_high = false
    elseif (a_mid - delt_a) .< minSpin
        run_low = false
    end

    if run_high
        a_list_high = LinRange(a_mid + delt_a, maxSpin, n_high)
        if !cheby
            fn = "rate_sve/Imag_erg_pos_$(state_str).npz"
        else
            fn = "rate_sve/Imag_ergC_pos_$(state_str).npz"
        end
        if isfile(fn)
            file_in = npzread(fn)
        
            file_in[:, :, 3] .*= 2.0
            alpha_load = unique(file_in[:,:,1])
            a_load = unique(file_in[:,:,2])
            file_in[file_in[:, :, 3] .<= 0.0, 3] .= 1e-100
            itp = LinearInterpolation((alpha_load, a_load), log10.(file_in[:, :, 3] ./ (GNew * M)), extrapolation_bc=Line())
            out_high = 10 .^itp(alph, a_list_high)
        else
            a_list_high = []
            out_high = []
            run_high = false
        end
    else
        a_list_high = []
        out_high = []
    end

    if run_low
        a_list_low = LinRange(minSpin, a_mid - delt_a, n_low)
        if !cheby
            fn = "rate_sve/Imag_erg_neg_$(state_str).npz"
        else
            fn = "rate_sve/Imag_ergC_neg_$(state_str).npz"
        end
        if isfile(fn)
            file_in = npzread(fn)
        
            file_in[:, :, 3] .*= 2.0
            alpha_load = unique(file_in[:,:,1])
            a_load = unique(file_in[:,:,2])
            file_in[file_in[:, :, 3] .<= 0.0, 3] .= 1e-100
            itp = LinearInterpolation((alpha_load, a_load), log10.(file_in[:, :, 3] ./ (GNew * M)), extrapolation_bc=Line())
            out_low = 10 .^itp(alph, a_list_low)
        else
            a_list_low = []
            out_low = []
            run_low = false
        end
    else
        a_list_low = []
        out_low = []
    end
    
    return a_mid, run_high, a_list_high, out_high, run_low, a_list_low, out_low

end

function check_slv_rad(alph, m)
    if (m == 1)&&(alph < 0.03)
        return true
    elseif (m == 2)&&(alph < 0.03)
        return true
    elseif (m==3)&&(alph < 0.16)
        return true
    elseif (m == 4)&&(alph < 0.40)
        return true
    else
        return false
    end
end

function precomputed_spin1(alph, a, M)
    rp = 1 .+ sqrt.(1 .- a.^2)
    Alist = [0.142 -1.17 3.024 -3.234 1.244; -1.298 6.601 -15.21 14.47 -5.070; 0.726 -8.516 15.43 -11.15 3.277]
    Blist = [-27.76 114.9 -311.1 177.2; -14.05 20.78 -36.84 58.37; 14.78 -4.574 -248.5 108.1]
    Clist = [48.86 -8.430 45.66 -132.8 52.48; -31.20 32.52 -73.50 161.0 -91.27; 189.1 -85.32 388.6 -1057 566.1]
    alphaL = []
    betaL = []

    for i in 1:3
        temp = 0.0
        for j in 1:length(Alist[i,:])
            temp += Alist[i,j] .* (1-a.^2).^((j-1)/2)
        end
        push!(alphaL, temp)
        
        temp = 0.0
        for j in 1:length(Clist[i,:])
            if j > 1
                temp += Blist[i,j-1] .* (1-a.^2).^((j-1)/2)
            end
            temp += Clist[i,j] .* a.^(j .- 1)
        end
        push!(betaL, temp)
    end
    
    wR = alph .* (1 .+ alphaL[1] .* alph .+ alphaL[2] .* alph.^2 .+ alphaL[3] .* alph.^3)
    wI = betaL[1] .*  alph .^ 7 .* (1 .+ betaL[2] .* alph .+ betaL[3] .* alph.^2 ) .* (a .- wR .* rp)
    
    return wR ./ (GNew .* M), wI ./ (GNew .* M)
end


function eigensys_Cheby(M, atilde, mu, n, l0, m; prec=200, L=4, Npoints=60, Iter=30, debug=false, return_wf=false, der_acc=1e-20, cvg_acc=1e-10, Npts_r=1000, nu_guess=nothing, return_nu=false, sfty_run=false)
    # L field spherical harmonic truncation l-eigenstate
    # Npoints number of Chebyshev interpolation points
    # Iter number of iterations for the non-linear inversion
    
    
    setprecision(BigFloat, prec)  # Higher than requested precision for calculations
    mu = BigFloat(mu)
    M = BigFloat(M)
    a = BigFloat(atilde) * M  # BH spin
    alph = BigFloat(mu .* GNew)
    n0 = n - l0 - 1 # field overtone number
   
    
    
    # Define auxiliary quantities
    rminus = M - sqrt(M^2 - a^2)
    rplus = M + sqrt(M^2 - a^2)

    # Helper function for Kronecker delta
    function kronecker_delta(i, j)
        i == j ? 1 : 0
    end

    # Defining ζ-r map
    function ζmap(r)
        return (r - sqrt(4*rplus*(r - rminus) + rminus^2))/(r - rminus)
    end

    function rmap(ζ)
        return (-rminus .+ 4 .* rplus .+ rminus .* ζ.^2) ./ (-1 .+ ζ).^2
    end

    # Black hole function
    function Δ(r)
        return r^2 - 2*M*r + a^2
    end


    # Define coupling function
    function c_coupling(l, j)
        if j-2 <= l <= j+2 && -l <= m <= l && -j <= m <= j
            return (1/3) * kronecker_delta(l, j) +
                   (2/3) * sqrt((2j + 1)/(2l + 1)) *
                   clebschgordan(j, m, 2, 0, l, m) *
                   clebschgordan(j, 0, 2, 0, l, 0)
        else
            return (1/3) * kronecker_delta(l, j)
        end
    end

    # Derivatives of the ζ-r map
    function Dζmap(r)
        return (rminus^2 + 2*r*rplus -
               rminus*(2*rplus + sqrt(rminus^2 + 4*r*rplus - 4*rminus*rplus)))/
               ((r - rminus)^2 * sqrt(rminus^2 + 4*r*rplus - 4*rminus*rplus))
    end

    function DDζmap(r)
        return 1/(r - rminus)^3 * 2 * (r +
               (2*(r - rminus)^2 * rplus^2)/(rminus^2 + 4*r*rplus - 4*rminus*rplus)^(3/2) -
               sqrt(rminus^2 + 4*r*rplus - 4*rminus*rplus) -
               (r - rminus) * (1 - (2*rplus)/sqrt(rminus^2 + 4*r*rplus - 4*rminus*rplus)))
    end

    # Define Chebyshev interpolation points
    ζ = [cos((BigFloat(π) * (2*k + 1))/(2 * (Npoints + 1))) for k in 0:Npoints]  # Chebyshev nodes
    w = [sin((2*BigFloat(π)*k*(Npoints + 2) + π)/(2 * (Npoints + 1))) for k in 0:Npoints]  # Chebyshev weights

    # Derivative matrices in polynomial base
    DerivP = zeros(BigFloat, Npoints+1, Npoints+1)
    for k in 0:Npoints, n in 0:Npoints
        if k != n
            DerivP[k+1, n+1] = w[k+1]/w[n+1]/(ζ[n+1] - ζ[k+1])
        else
            sum_val = sum(w[kk+1]/w[n+1]/(ζ[n+1] - ζ[kk+1]) for kk in 0:Npoints if kk != n)
            DerivP[k+1, n+1] = -sum_val
        end
    end

    # Second derivative matrix
    DDerivP = zeros(BigFloat, Npoints+1, Npoints+1)
    for k in 0:Npoints, n in 0:Npoints
        if k != n
            DDerivP[k+1, n+1] = 2 * DerivP[k+1, n+1] * DerivP[n+1, n+1] -
                               (2 * DerivP[k+1, n+1])/(ζ[n+1] - ζ[k+1])
        else
            sum_val = sum(2 * DerivP[kk+1, n+1]/(ζ[n+1] - ζ[kk+1]) for kk in 0:Npoints if kk != n)
            DDerivP[k+1, n+1] = 2 * DerivP[k+1, n+1] * DerivP[n+1, n+1] + sum_val
        end
    end

    function generalized_laguerre(n, α, x)
        return (-1).^ n ./ factorial(big(n)) .* HypergeometricFunctions.U.(-n, α+1 ,x)
        # return binomial(n + α, n) .* HypergeometricFunctions.M.(-n, α+1 , x)
    end


    # Function for hydrogenic frequency parameter
    function calc_Nν_initial()
        return l0 + 1 + n0 +
            im * (alph * M)^(4*l0 + 2) * ((a*m)/M - 2*alph*rplus) *
            (2^(4*l0 + 2) * factorial(big(2*l0 + 1 + n0)))/
            ((l0 + 1 + n0)^(2*l0 + 1) * factorial(big(n0))) *
            (factorial(big(l0))/(factorial(big(2*l0)) * factorial(big(2*l0 + 1))))^2 *
            prod(j^2 * (1 - a^2/M^2) + ((a*m)/M - 2*alph*rplus)^2 for j in 1:l0)
    end
    


    # Now we can define ω
    function calc_ω(ν_val)
        return alph .* sqrt.(1 .- (alph .^2 .* M.^2) ./ ν_val.^2)
    end
    function calc_ω_inv(ω)
        return sqrt.(((1 - (ω ./ alph).^2 ) ./ (alph^2*M^2) ).^-1)
    end
    
    # Hydrogenic frequency parameter ν (initial value)
    if isnothing(nu_guess)
        Nν = calc_Nν_initial()
    else
        erg_in = BigFloat(real(nu_guess)) .+ im * BigFloat(imag(nu_guess))
        Nν = calc_ω_inv(erg_in)
        
    end
    
    function R(r, ν_val)
        # Hydrogenic radial wave function
        return ((2*r*M*alph^2)/(l0 + n0 + 1))^l0 *
                exp(-((r*M*alph^2)/(l0 + n0 + 1))) *
                generalized_laguerre(n - l0 - 1, 2 * l0 + 1, (2*r*M*alph^2)/(l0 + n0 + 1))
    end
    
    if debug
        println("starting erg \t", calc_ω(Nν).* M)
    end

    # Define P and A values as functions of ν
    function calc_Pplus(ν_val)
        ω_val = calc_ω(ν_val)
        return (m*a - 2*ω_val*M*rplus)/(rplus - rminus)
    end

    function calc_Pminus(ν_val)
        ω_val = calc_ω(ν_val)
        return (m*a - 2*ω_val*M*rminus)/(rplus - rminus)
    end

    function calc_Aplus(ν_val)
        ω_val = calc_ω(ν_val)
        return 2*((m*a - 2*ω_val*M^2)/(rplus - rminus))^2 -
               4*ω_val^2*M^2 + (alph^2 - 2*ω_val^2)*M*(rplus - rminus) +
               (2*M^2 - a^2)*(alph^2 - ω_val^2)
    end

    function calc_Aminus(ν_val)
        ω_val = calc_ω(ν_val)
        return 2*((m*a - 2*ω_val*M^2)/(rplus - rminus))^2 -
               4*ω_val^2*M^2 - (alph^2 - 2*ω_val^2)*M*(rplus - rminus) +
               (2*M^2 - a^2)*(alph^2 - ω_val^2)
    end

    # Define asymptotic behavior functions as functions of ν
    function F(r, ν_val)
        ω_val = calc_ω(ν_val)
        Pplus_val = calc_Pplus(ν_val)
        return ((r - rplus)/(r - rminus))^(im * Pplus_val) *
               (r - rminus)^(-1 + ν_val - (2*alph^2*M^2)/ν_val) *
               exp(-(alph^2*M)/ν_val * (r - rplus))
    end

    function DF(r, ν_val)
        ω_val = calc_ω(ν_val)
        Pplus_val = calc_Pplus(ν_val)
        return -(1/ν_val) * exp(-((M*(r - rplus)*alph^2)/ν_val)) *
               (r - rminus)^(-3 - (2*M^2*alph^2)/ν_val + ν_val) *
               ((r - rplus)/(r - rminus))^(-1 + im*Pplus_val) *
               (2*M^2*(r - rplus)*alph^2 + M*(r - rminus)*(r - rplus)*alph^2 +
                ν_val*(r + im*Pplus_val*(rminus - rplus) + rplus*(-1 + ν_val) - r*ν_val))
    end

    function DDF(r, ν_val)
        ω_val = calc_ω(ν_val)
        Pplus_val = calc_Pplus(ν_val)
        return (1/((r - rplus)^2 * ν_val^2)) * exp(-((M*(r - rplus)*alph^2)/ν_val)) *
               (r - rminus)^(-3 - (2*M^2*alph^2)/ν_val + ν_val) *
               ((r - rplus)/(r - rminus))^(im*Pplus_val) *
               (4*M^4*(r - rplus)^2*alph^4 + 4*M^3*(r - rminus)*(r - rplus)^2*alph^4 -
                2*M*(r - rminus)*(r - rplus)*alph^2*ν_val*(-im*Pplus_val*(rminus - rplus) +
                rplus + r*(-1 + ν_val) - rplus*ν_val) + M^2*(r - rplus)*alph^2*(r^3*alph^2 -
                rminus^2*rplus*alph^2 - r^2*(2*rminus + rplus)*alph^2 + 4*im*Pplus_val*rminus*ν_val +
                2*rplus*ν_val*(-3 - 2*im*Pplus_val + 2*ν_val) + r*(rminus^2*alph^2 + 2*rminus*rplus*alph^2 +
                2*(3 - 2*ν_val)*ν_val)) + ν_val^2*(-Pplus_val^2*(rminus - rplus)^2 -
                im*Pplus_val*(rminus - rplus)*(rminus + 3*rplus - 2*rplus*ν_val) +
                2*r*(-2 + ν_val)*(-im*Pplus_val*(rminus - rplus) + rplus - rplus*ν_val) +
                r^2*(2 - 3*ν_val + ν_val^2) + rplus^2*(2 - 3*ν_val + ν_val^2)))
    end



    # Function to build equation coefficients
    function build_coefficients(ν_val)
        ω_val = calc_ω(ν_val)
        Pplus_val = calc_Pplus(ν_val)
        Pminus_val = calc_Pminus(ν_val)
        Aplus_val = calc_Aplus(ν_val)
        Aminus_val = calc_Aminus(ν_val)
        
        C1 = zeros(Complex{BigFloat}, L+1, Npoints+1)
        C2 = zeros(Complex{BigFloat}, L+1, Npoints+1)
        C3plus = zeros(Complex{BigFloat}, L+1, Npoints+1)
        C3minus = zeros(Complex{BigFloat}, L+1, Npoints+1)
        
        
        gam = im * a * sqrt.(ω_val.^2 .- alph.^2)
        
        for l in 0:L, n in 0:Npoints
            LLM = l * (l + 1)
            LLM += (-1 + 2 * l * (l + 1) - 2 * m.^2) * gam.^2 ./ (-3 + 4 * l * (l + 1))
            LLM += ((l - m - 1 * (l - m) * (l + m) * (l + m - 1)) ./ ((-3 + 2 * l) * (2 * l - 1).^2) - (l + 1 - m) * (2 * l - m) * (l + m + 1) * (2 + l + m) ./ ((3 + 2 * l).^2 * (5 + 2 * l))) * gam.^4 ./ (2 * (1 + 2 * l))
            LLM += (4 * ((-1 + 4 * m^2) * (l * (1 + l) * (121 + l * (1 + l) * (213 + 8 * l * (1 + l) * (-37 + 10 * l * (1 + l)))) - 2 * l * (1 + l) * (-137 + 56 * l * (1 + l) * (3 + 2 * l * (1 + l))) * m^2 + (705 + 8 * l * (1 + l) * (125 + 18 * l * (1 + l))) * m^4 - 15 * (1 + 46 * m^2))) * gam^6) / ((-5 + 2 * l) * (-3 + 2 * l) * (5 + 2 * l) * (7 + 2 * l) * (-3 + 4 * l * (1 + l))^5)
            
            r_n = rmap(ζ[n+1])
            C1[l+1, n+1] = (1/(r_n - rplus) + 1/(r_n - rminus)) * 1/Dζmap(r_n) +
                           (2*DF(r_n, ν_val))/F(r_n, ν_val) * 1/Dζmap(r_n) +
                           DDζmap(r_n)/(Dζmap(r_n)^2)
            
            C2[l+1, n+1] = 1/(Dζmap(r_n))^2 * DDF(r_n, ν_val)/F(r_n, ν_val) +
                           (1/(r_n - rplus) + 1/(r_n - rminus)) *
                           1/(Dζmap(r_n))^2 * DF(r_n, ν_val)/F(r_n, ν_val) +
                           1/(Dζmap(r_n))^2 * (Pplus_val^2/(r_n - rplus)^2 +
                           Pminus_val^2/(r_n - rminus)^2 -
                           Aplus_val/((rplus - rminus)*(r_n - rplus)) +
                           Aminus_val/((rplus - rminus)*(r_n - rminus)) - (alph^2 - ω_val^2)) -
                           1/(Dζmap(r_n))^2 * 1/Δ(r_n) * LLM
                           # 1/(Dζmap(r_n))^2 * 1/Δ(r_n) .* (l*(l + 1) + a^2*c_coupling(l, l)*(alph^2 - ω_val^2))
                           
            
            C3plus[l+1, n+1] = -(1/(Dζmap(r_n))^2) * 1/Δ(r_n) *
                               a^2 * (alph^2 - ω_val^2) * c_coupling(l, l + 2)
            
            C3minus[l+1, n+1] = -(1/(Dζmap(r_n))^2) * 1/Δ(r_n) *
                                a^2 * (alph^2 - ω_val^2) * c_coupling(l, l - 2)
        end
        
        return C1, C2, C3plus, C3minus
    end

    # Function to build the system matrix
    function build_matrix(ν_val)
        C1, C2, C3plus, C3minus = build_coefficients(ν_val)
        
        # Size of the full matrix
        N = (L+1)*(Npoints+1)
        MatrixOp = zeros(Complex{BigFloat}, N, N)
        
        for l in 0:L
            for n in 0:Npoints
                row = l*(Npoints+1) + n + 1
                
                # Diagonal part (l,l)
                for k in 0:Npoints
                    col = l*(Npoints+1) + k + 1
                    MatrixOp[row, col] = DDerivP[k+1, n+1] + C1[l+1, n+1]*DerivP[k+1, n+1]
                    if k == n
                        MatrixOp[row, col] += C2[l+1, n+1]
                    end
                end
                
                # Off-diagonal part (l,l-2)
                if l >= 2
                    col = (l-2)*(Npoints+1) + n + 1
                    MatrixOp[row, col] = C3minus[l+1, n+1]
                end
                
                # Off-diagonal part (l,l+2)
                if l+2 <= L
                    col = (l+2)*(Npoints+1) + n + 1
                    MatrixOp[row, col] = C3plus[l+1, n+1]
                end
            end
        end
        
        return MatrixOp
    end
    
    
    
    # Function to build the derivative of the system matrix with respect to ν
    function build_derivative_matrix(ν_val)
        # This is a complex operation, would need symbolic differentiation
        # Placeholder implementation - actual implementation would need careful derivative calculations
        
                
        
        # Use numerical differentiation as an approximation
        ε = der_acc
        del_nu = ν_val .* ε
        """
        # order x^2 NOT HIGH ENOUGH
        # testD = (build_matrix(ν_val .* (1 .+ ε)) - build_matrix(ν_val .* (1 .- ε)))/(2*ν_val*ε)
        # x^4
        der_out =  (-build_matrix(ν_val .+ 2 .* del_nu) .+ 8 .* build_matrix(ν_val .+ del_nu)  .- 8 .* build_matrix(ν_val .- del_nu) .+ build_matrix(ν_val .- 2 .* del_nu))/(12 .* del_nu)
        
        
        # x^6
        der_out = -(  build_matrix(ν_val .- 3 .* del_nu)
               .- 9 .* build_matrix(ν_val .- 2 .* del_nu)
               .+ 45 .* build_matrix(ν_val .- 1 .* del_nu)
               .- 45 .* build_matrix(ν_val .+ 1 .* del_nu)
               .+ 9 .* build_matrix(ν_val .+ 2 .* del_nu)
               .- 1 .* build_matrix(ν_val .+ 3 .* del_nu)
              ) ./ (60 .* del_nu)
              
        # x^8
        der_out = ( -3  .* build_matrix(ν_val .+ 4 .* del_nu)
               +32  .* build_matrix(ν_val .+ 3 .* del_nu)
               -168 .* build_matrix(ν_val .+ 2 .* del_nu)
               +672 .* build_matrix(ν_val .+ 1 .* del_nu)
               -672 .* build_matrix(ν_val .- 1 .* del_nu)
               +168 .* build_matrix(ν_val .- 2 .* del_nu)
               -32  .* build_matrix(ν_val .- 3 .* del_nu)
               +3   .* build_matrix(ν_val .- 4 .* del_nu)
              ) ./ (840 .* del_nu)
        """
        # der_out =  (-build_matrix(ν_val .+ 2 .* del_nu) .+ 8 .* build_matrix(ν_val .+ del_nu)  .- 8 .* build_matrix(ν_val .- del_nu) .+ build_matrix(ν_val .- 2 .* del_nu))/(12 .* del_nu)
        
        der_out = ( -3  .* build_matrix(ν_val .+ 4 .* del_nu)
               +32  .* build_matrix(ν_val .+ 3 .* del_nu)
               -168 .* build_matrix(ν_val .+ 2 .* del_nu)
               +672 .* build_matrix(ν_val .+ 1 .* del_nu)
               -672 .* build_matrix(ν_val .- 1 .* del_nu)
               +168 .* build_matrix(ν_val .- 2 .* del_nu)
               -32  .* build_matrix(ν_val .- 3 .* del_nu)
               +3   .* build_matrix(ν_val .- 4 .* del_nu)
              ) ./ (840 .* del_nu)

    
        return der_out
    end

    # Initialize solution vectors for non-linear inverse iteration
    Nν_values = Vector{Complex{BigFloat}}(undef, Iter+1)
    Nν_values[1] = Nν  # Initial value

    # Function to build initial guess vector
    function build_initial_vector()
        N = (L+1)*(Npoints+1)
        v = zeros(Complex{BigFloat}, N)
        
        counter = 1
        for j in 0:L
            for k in 0:Npoints
                v[counter] = (kronecker_delta(j, l0) * R(rmap(ζ[k+1]), Nν))/F(rmap(ζ[k+1]), Nν)
                counter += 1
            end
        end
        
        return v
    end

    # Normalization vector for iteration
    function build_normalization_vector()
        N = (L+1)*(Npoints+1)
        return [exp(im * jk * π/N)/sqrt(N) for jk in 1:N]
    end

    # Perform non-linear inverse iteration
    function run_nonlinear_inverse_iteration()
        # Initialize vectors
        Bnum = Vector{Vector{Complex{BigFloat}}}(undef, Iter+1)
        BnumNorm = Vector{Vector{Complex{BigFloat}}}(undef, Iter+1)
        
        # Initial guess
        
        Bnum[1] = build_initial_vector()
        
        BnumNorm[1] = Bnum[1] / sqrt(dot(conj(Bnum[1]), Bnum[1]))
        
        # Normalization vector
        VEC = build_normalization_vector()
        
        # Iteration
        final_idx = 1
        i = 1
        while i < (Iter + 1)

            A_i = build_matrix(Nν_values[i])
            b_i = build_derivative_matrix(Nν_values[i]) * BnumNorm[i]
            
            # Solve linear system
            Bnum[i+1] = A_i \ b_i
            
            # Update eigenvalue estimate
            Nν_values[i+1] = Nν_values[i] - dot(conj(VEC), BnumNorm[i]) / dot(conj(VEC), Bnum[i+1])
            

            # Normalize solution vector
            BnumNorm[i+1] = Bnum[i+1] / sqrt(dot(conj(Bnum[i+1]), Bnum[i+1]))
            
            if debug
                println("Iteration $i: ν = $(Nν_values[i+1])")
                if isnan(Nν_values[i+1])
                    println("HERE \n")
                    prec += 50
                    setprecision(BigFloat, prec)
                end
             end
            
            if i > 1
                erg_test_init = calc_ω(Nν_values[i])
                erg_test = calc_ω(Nν_values[i+1])
                test_convR = abs.(real.(erg_test .- erg_test_init) ./ minimum([real.(erg_test_init), real.(erg_test)]))
                test_convI = abs.(imag.(erg_test .- erg_test_init) ./ minimum([imag.(erg_test_init), imag.(erg_test)]))
                if (test_convR < cvg_acc)&&(test_convI < cvg_acc)
                    final_idx = i
                    i = Iter + 1
                end
            end
            i += 1
            
        end
        
        return Nν_values, Bnum, BnumNorm, final_idx
    end


    # Run the algorithm
    Nν_values, Bnum, BnumNorm, final_idx = run_nonlinear_inverse_iteration()
    erg_out = calc_ω(Nν_values[final_idx]) .* M
    

    
    if debug
        println("Final ν value: ", erg_out)
    end
    
    if sfty_run
        sfty = false
        Npoints += 15
        while !sfty
            eR, eI = eigensys_Cheby(M, atilde, mu, n, l0, m; prec=prec, L=L, Npoints=Npoints, Iter=Iter, debug=debug, return_wf=false, der_acc=der_acc, cvg_acc=cvg_acc, Npts_r=Npts_r, nu_guess=(erg_out ./ M), return_nu=false, sfty_run=false)
            imdiff = abs.(log10.(abs.(eI)) .- log10.(abs.(imag(erg_out)))) ./ minimum([abs.(log10.(abs.(eI))),abs.(log10.(abs.(imag(erg_out))))])
            realdiff = abs.(log10.(abs.(eR)) .- log10.(abs.(real(erg_out)))) ./ minimum([abs.(log10.(abs.(eR))),abs.(log10.(abs.(real(erg_out))))])
            if debug
                println(Npoints, "\t", eI, "\t", imag(erg_out))
            end
            if (imdiff .<= 0.1)&&(realdiff .<= 1e-3)
                erg_out = eR + im * eI
                sfty = true
            else
                erg_out = eR + im * eI
            end
            Npoints += 15
            if Npoints > 150
                println("Npoints exceeded 150.... ")
                erg_out = eR + im * eI
                sfty = true
            end
        end
    end
    
    
    if !return_wf
        if !return_nu
            return real(erg_out), imag(erg_out) # real(G * M * omega), imag(G * M * omega)
        else
            return real(erg_out), imag(erg_out), Nν_values[final_idx] # real(G * M * omega), imag(G * M * omega),
        end
    end
    
    # Calculate ω based on the final Nν value from iteration
    ω = alph * sqrt(1 - (alph^2 * M^2)/(Nν_values[final_idx])^2)
    kk = sqrt.(ω.^2 .- alph.^2)
    
    
    # Calculate Pplus for the current ν value
    Pplus = calc_Pplus(Nν_values[final_idx])
    Pmns = calc_Pminus(Nν_values[final_idx])
    beta = 2 * im .* Pplus
    gam = -2 .* im .* Pmns
    alpha_other = 2 * im * kk .* (rplus .- rminus)
    deltt = calc_Aplus(Nν_values[final_idx]) .- calc_Aminus(Nν_values[final_idx])
    gam_llm = im * a * sqrt.(ω.^2 .- alph.^2)
    LLM = l0 * (l0 + 1)
    LLM += (-1 + 2 * l0 * (l0 + 1) - 2 * m.^2) * gam_llm.^2 ./ (-3 + 4 * l0 * (l0 + 1))
    LLM += ((l0 - m - 1 * (l0 - m) * (l0 + m) * (l0 + m - 1)) ./ ((-3 + 2 * l0) * (2 * l0 - 1).^2) - (l0 + 1 - m) * (2 * l0 - m) * (l0 + m + 1) * (2 + l0 + m) ./ ((3 + 2 * l0).^2 * (5 + 2 * l0))) * gam_llm.^4 ./ (2 * (1 + 2 * l0))
    LLM += (4 * ((-1 + 4 * m^2) * (l0 * (1 + l0) * (121 + l0 * (1 + l0) * (213 + 8 * l0 * (1 + l0) * (-37 + 10 * l0 * (1 + l0)))) - 2 * l0 * (1 + l0) * (-137 + 56 * l0 * (1 + l0) * (3 + 2 * l0 * (1 + l0))) * m^2 + (705 + 8 * l0 * (1 + l0) * (125 + 18 * l0 * (1 + l0))) * m^4 - 15 * (1 + 46 * m^2))) * gam_llm^6) / ((-5 + 2 * l0) * (-3 + 2 * l0) * (5 + 2 * l0) * (7 + 2 * l0) * (-3 + 4 * l0 * (1 + l0))^5)
    eta = - (Pplus.^2 .+ Pmns.^2 .+ LLM .+ alph.^2 .* rplus.^2 .- ω.^2 .* (4 .* M.^2 .+ 2 .* M .* rplus .+ rplus.^2))
    phi = (alpha_other .- beta .- gam .+ alpha_other .* beta .- beta .* gam) ./ 2 .- eta
    nuN = (alpha_other .+ beta .+ gam .+ alpha_other .* gam  .+ beta .* gam) ./ 2 .+ eta .+ deltt
    
    ### choose r values i want...
    r_max_shrt = 2.0 .^(2.0 .* n .- 2 .* (1 .+ n)) .* gamma(2 .+ 2 .* n) ./ alph.^2 ./ factorial(big(2 .* n - 1)) .* 3.0
    r_max = n ./ alph.^2 .* 100.0
    if (R(r_max, Nν_values[final_idx]) ./ R(rplus  .* (1.0 .+ 1e-3), Nν_values[final_idx])) .> 1
        r_max *= 10
    end
    r_vals = rlist = 10 .^(range(log10.(rplus  .* (1.0 .+ 1e-3)), log10.(r_max), Npts_r))
    rout_temp = R.(r_vals, Nν_values[final_idx])
    rout_temp = abs.(rout_temp)
    rout_temp ./= trapz(rout_temp.^2 .* r_vals.^2, r_vals)
    r_thresh = nothing
    for i in 1:length(r_vals)
        if r_vals[end - i + 1] .< r_max_shrt
            r_thresh = r_vals[end - i + 1]
            break
        end
    end

    zz = -(r_vals .- rplus) ./ (rplus .- rminus)
    phinu = Float64(real.(phi + nuN)) .+ im .* Float64(imag.(phi + nuN))
    
    zz_short = zz[abs.(zz) .< r_thresh] # cant call large zz values of heun (takes forever...)
    heunc = solve_heun_switch(phi, phi + nuN, 1.0 + beta, 1.0 + gam, alpha_other, zz_short)

    
    y_values_short = exp.(-im .* kk .* (r_vals[abs.(zz) .< r_thresh] .- rplus)) .* zz_short.^(im .* Pplus) .* (zz_short .- 1.0).^(- im .* Pmns) .* heunc
    
    y_values = []
    rmax_cut = r_vals[abs.(zz) .< r_thresh][end]
    
    # work backward to find where the solution fails
    imax = length(zz_short)
    if abs.(y_values_short[end]) .> abs.(y_values_short[end-1])
        for i in 1:(length(zz_short)-1)
            if abs.(y_values_short[end-i+1]) .> abs.(y_values_short[end-i])
                y_values_short[end-i+1] = 0.0 + im * 0.0
            else
                imax = length(zz_short) - i - 1
                break
            end
        end
    end
    # now go a little further back because you didn't go far enough
    imax = argmin(abs.(r_vals .- r_vals[imax] .* 0.9))
    r_max_new = r_vals[imax]
    
    for i in 1:Npts_r
        if r_vals[i] <= r_max_new
            push!(y_values, y_values_short[i])
        else
            push!(y_values, y_values_short[imax] .* exp.(im .* kk .* (r_vals[i] .- r_max_new)) .* (r_vals[i] ./ r_max_new).^2 )
        end
    end
    rlist = r_vals ./ M
    
    
    nm2 = trapz(y_values .* conj.(y_values) .* rlist.^2, rlist)
    y_values ./= sqrt.(nm2)
    
    number_max_cnt = count_local_maxima(real.(y_values .* conj.(y_values)))
    if number_max_cnt != (n - l0)
        println("Nmax count :", number_max_cnt)
        y_values = rout_temp
        nm2 = trapz(y_values .* conj.(y_values) .* rlist.^2, rlist)
        y_values ./= sqrt.(nm2)
    end
    
    if !return_nu
        return real(erg_out), imag(erg_out), rlist, y_values # everything normalized by GM, radial WF may not resolve at large r
    else
        return real(erg_out), imag(erg_out), rlist, y_values, Nν_values[final_idx]
    end
end

function count_local_maxima(arr::AbstractVector)
    count = 0
    for i in 2:length(arr)-1
        if arr[i] > arr[i-1] && arr[i] > arr[i+1]
            count += 1
        end
    end
    return count
end

function trapz(y,x)
    return sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
end

function normalize_mathlink_str(s::String)
    s = strip(s)
    # If Mathematica prints something like 1.23`15.7*^309 or 1.23`15.7*^309
    # remove the backtick precision part (e.g. `15.7...) and convert *^ -> e
    s = replace(s, r"`[^*e]*" => "")   # remove backtick + subsequent chars until * or e
    s = replace(s, "*^" => "e")
    return s
end


# Strict reconstruction if you want to avoid parse of exponent form:
function reconstruct_from_parts(mantissa_str::AbstractString, exponent::Integer)
    ndigits = count(isdigit, mantissa_str)
    bits = ceil(Int, ndigits * 3.321928094887362) + 64
    setprecision(bits) do
        mant = parse(BigFloat, mantissa_str)
        return mant * big(10)^exponent
    end
end

# ---------------------------
# Top-level converter for MathLink.WReal (or similar)
# ---------------------------

"""
    mlwreal_to_bigfloat(w)

Convert a MathLink.WReal-like object `w` to a `BigFloat`.

Strategy:
- try `string(w)` / `repr(w)` and parse the returned representation;
- if `w` exposes `mantissa` and `exponent` properties, use those to reconstruct exactly;
- otherwise return an informative error with a dump of the object.
"""
function mlwreal_to_bigfloat(w)
    # 1) try string representation (most wrappers provide a readable string)
    s = nothing
    try
        s = string(w)
    catch err
        try
            s = repr(w)
        catch err
            s = nothing
        end
    end

    if s !== nothing && (occursin("*^", s) || occursin("e", s) || occursin("`", s))
        # likely a Mathematica numeric form; normalize & parse
        return parse_mathlink_number(s)
    end

    # 2) Some wrappers return a small struct with mantissa/exponent fields
    #    Try to use them if present.
    if hasproperty(w, :mantissa) && hasproperty(w, :exponent)
        mant = getproperty(w, :mantissa)
        exp  = getproperty(w, :exponent)
        # mantissa may itself be numeric or string; coerce to string
        mantstr = isa(mant, AbstractString) ? mant : string(mant)
        return reconstruct_from_parts(mantstr, Int(exp))
    end

    # 3) Try other common names
    for name in (:mant, :mantissa_str, :man, :m, :exp, :exponent_int, :e)
        if hasproperty(w, name)
            try
                mant = getproperty(w, name)
                exp = hasproperty(w, :exponent) ? getproperty(w, :exponent) :
                      (hasproperty(w, :exp) ? getproperty(w, :exp) : nothing)
                if exp !== nothing
                    mantstr = isa(mant, AbstractString) ? mant : string(mant)
                    return reconstruct_from_parts(mantstr, Int(exp))
                end
            catch
                # ignore and continue
            end
        end
    end

    # 4) last resort: if string is something like "WReal( ... )", try to extract inner "*^" substring
    if s !== nothing
        m = match(r"([+-]?\d*\.?\d+(?:`[^*]*)?\*\^[+-]?\d+)", s)
        if m !== nothing
            return parse_mathlink_number(m.captures[1])
        end
    end

    error("Unable to convert MathLink object to BigFloat. Representation: $(s === nothing ? repr(w) : s)")
end

function parse_mathlink_number(w)
    # Step 1: get a string representation
    s = string(w)

    # Step 2: remove Mathematica-style precision and exponent notation
    s = strip(s)
    s = replace(s, r"`[^*e]*" => "")   # remove backtick precision (e.g. `15.7...)
    s = replace(s, "*^" => "e")        # replace Mathematica exponent with 'e'

    # Step 3: choose precision based on digits
    m = match(r"^([+-]?\d*\.?\d+)", s)
    mantissa = m === nothing ? "1.0" : m.captures[1]
    ndigits = count(isdigit, mantissa)
    bits = ceil(Int, ndigits * 3.321928) + 64

    # Step 4: parse into BigFloat
    setprecision(bits) do
        return parse(BigFloat, s)
    end
end

# Gaunt coefficient G(l1,l2,l3; m1,m2,m3)
function gaunt(l1::Integer, l2::Integer, l3::Integer,
               m1::Integer, m2::Integer, m3::Integer)
    # quick zero checks
    if m1 + m2 + m3 != 0
        return 0.0
    end
    # triangle condition for 3j with ms = 0
    if !((abs(l1-l2) <= l3 <= l1 + l2) && ((l1 + l2 + l3) % 1 == 0))
        return 0.0
    end
    if (m1 > l1)||(m2 > l2)||(m3>l3)
        return 0.0
    end
    pref = sqrt((2l1+1)*(2l2+1)*(2l3+1) / (4 * pi))
    t1 = Float64.(wigner3j(l1, l2, l3, 0, 0, 0))
    t2 = Float64.(wigner3j(l1, l2, l3, m1, m2, m3))
    return pref * t1 * t2
end

"""
    integral4(l1,m1,l2,m2,l3,m3,l4,m4; tol=1e-12)

Compute I = ∫ Y_{l1 m1} Y_{l2 m2} Y_{l3 m3} Y_{l4 m4} dΩ
Returns (value, is_zero::Bool) where is_zero tests abs(value) < tol.
"""
function integral4(l1::Integer, m1::Integer,
                   l2::Integer, m2::Integer,
                   l3::Integer, m3::Integer,
                   l4::Integer, m4::Integer; tol=1e-12)

    # Selection rule: total m must be zero
    if m1 + m2 + m3 + m4 != 0
        return 0.0, true
    end

    total = 0.0
    # L must satisfy triangle with (l1,l2) and with (l3,l4)
    Lmin = max(abs(l1 - l2), abs(l3 - l4))
    Lmax = min(l1 + l2, l3 + l4)
    for L in Lmin:Lmax
        # For a given L, only M = -(m1+m2) contributes (because Gaunt nonzero only when m1+m2+M=0, and similarly)
        M = -(m1 + m2)
        if abs(M) > L
            continue
        end
        g1 = gaunt(l1, l2, L, m1, m2, M)
        if g1 == 0.0
            continue
        end
        g2 = gaunt(l3, l4, L, m3, m4, -M)
        if g2 == 0.0
            continue
        end
        total += g1 * g2
    end

    is_zero = abs(total) < tol
    return total, is_zero
end

# test run of system
# @time wR, wI, rl, r3 = eigensys_Cheby(1, 0.9, 0.8557577459790169 ./ GNew, 4, 3, 3, debug=true, return_wf=true, L=4, Npoints = 70, Iter = 20,  der_acc=1e-20, cvg_acc=1e-10, prec=200, Npts_r=2000, sfty_run=false)

# @time wR, wI, rl, r3 = eigensys_Cheby(1, 0.9, 0.03 ./ GNew, 2, 1, 1, debug=true, return_wf=true, L=4, Npoints = 70, Iter = 20,  der_acc=1e-20, cvg_acc=1e-10, prec=200, Npts_r=200, sfty_run=false)

# @time wR, wI = find_im_part(0.3 ./ GNew, 1, 0.99, 2, 1, 1; debug=true, return_both=true, for_s_rates=true, Ntot_force=20000)
# println(wR ./ 0.3, "\t", wI ./ 0.3)
# println(wI)
# writedlm("test.dat", cat(real.(rl), real.(abs.(r3 .* conj(r3))), dims=2))

# M=1
# a=0.95
# alp = 0.05
# mu = alp ./ (M * GNew)
# @time wR, wI = find_im_part(mu, M, a, 3, 2, 0; debug=true, iter=500, xtol=1e-30, ftol=1e-90, return_both=true, for_s_rates=true, QNM=false, Ntot_force=10000)
#println(wR, "\t", wI)
# rl, r1, erg_1 = solve_radial(mu, M, a, 5, 4, 4; rpts=2000, return_erg=true, Ntot_safe=10000, use_heunc=true)
# writedlm("test_store/test_chb.dat", cat(real.(rl), real.(abs.(r1 .* conj(r1))), dims=2))

# wR, wI, rl, r3 = eigensys_Cheby(1, a, mu, 5, 4, 4, debug=true, return_wf=true, L=4, Npoints = 90, Iter = 20,  der_acc=1e-20, cvg_acc=1e-10, prec=200, Npts_r=200, sfty_run=false)
# writedlm("test_store/test_ctrue.dat", cat(real.(rl), real.(abs.(r3 .* conj.(r3))), dims=2))
