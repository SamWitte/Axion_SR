using SpecialFunctions
using Random
using OrdinaryDiffEq
using Statistics
using Distributions
using ForwardDiff: gradient, derivative, Dual, Partials, hessian
using DelimitedFiles
using NLsolve
using DifferentialEquations
using Interpolations
using HypergeometricFunctions
using SpinWeightedSpheroidalHarmonics
# using MCIntegration
# using Dates
# using LinearAlgebra
using Optim
include("Constants.jl")


function ergL(n, l, m, massB, MBH, a; full=true)
    # Key to next level
    alph = GNew * MBH * massB
    if full
        if l > 0
            return massB .* (1.0 .- alph.^2 ./ (2 .* n.^2) - alph.^4 ./ (8 * n.^4) + alph.^4 ./ n^4 .* (2 * l - 3 * n + 1) ./ (l + 0.5) + 2 * a * m * alph.^5 ./ n^3 ./ (l * (l + 0.5) * (l+1)))
        else
            return massB .* (1.0 .- alph.^2 ./ (2 .* n.^2) - alph.^4 ./ (8 * n.^4) + alph.^4 ./ n^4 .* (2 * l - 3 * n + 1) ./ (l + 0.5))
        end
    else
        return massB .* (1.0 .- alph.^2 ./ (2 .* n.^2) - alph.^4 ./ (8 * n.^4))
        # return massB .* (1.0 .- alph.^2 ./ (2 .* n.^2))
    end
end


function sr_rates(n, l, m, massB, MBH, aBH; impose_low_cut=0.001, solve_322=true)
    if (n==3)&&(l==2)&&(m==1)&&(solve_322==false)
        return 0.0
    end
    

    alph = GNew .* MBH .* massB
    
#    if (alph ./ l < impose_low_cut)&&(MBH < 1e2)
#        # We expect binaries to disrupt.
#        return 0.0
#    end
    
    
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
    Anl *= (factorial(Int(l)) ./ (factorial(Int(2 .* l)) .* factorial(Int(2 .* l .+ 1)))).^2
    Chilm = 1.0
    # erg = ergL(n, l, m, massB, MBH, aBH)
    for k in 1:Int(l)
        Chilm *= (k.^2 .* (1.0 .- aBH.^2) .+ (aBH * m .- 2 .* (rP ./ (GNew .* MBH)) .* alph).^2)
        # Chilm *= (k.^2 .* (1.0 .- aBH.^2) .+ (aBH * m .- 2 .* (rP ./ (GNew .* MBH)) .* alph).^2)
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
        rF[i] = 2 .* k .* exp.(pi ./ (2 .* k .* a0)) .* abs.(gamma(l + 1 - im ./ (k .* a0))) ./ factorial(2 * l + 1) .* (2 .* k .* r[i]).^l .* exp.(-im * k * r[i]) .* pFq((im ./ (k .* a0) + l + 1,), (2 * l + 2, ), 2 * im .* k .* r[i])
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
    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
    
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
    
    
    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
        

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
            if erg_4G != erg_ind
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
    

    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
    
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

function solve_radial(mu, M, a, n, l, m; rpts=1000, rmaxT=50, debug=false, iter=500, xtol=1e-20, ftol=1e-90, sve=false, fnm="test_store/WF_", return_erg=false, Ntot_safe=5000, eps_r=1e-10, QNM=false, QNM_ergs=nothing)
    ### dolan 2007
    # everything normalized
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
        wR, wI = find_im_part(mu, M, a, n, l, m; debug=debug, iter=iter, xtol=xtol, ftol=ftol, return_both=true, for_s_rates=true, QNM=false, Ntot_force=Ntot_safe)
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
        return rlist, zeros(length(rlist))
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
    
    # Normalize integration
    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
    
    
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
    
    
    
    
    
    # print(Float64.(real(Rout ./ sqrt.(nm))), "\n")
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
    
    # normalize solution!
    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
  
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
    Zlm = spin_weighted_spheroidal_harmonic(0, l, m, a .* erg)
    return Zlm
end

function find_im_part(mu, M, a, n, l, m; debug=false, Ntot_force=200, iter=10000, xtol=1e-20, ftol=1e-90, return_both=false, for_s_rates=false, QNM=false, QNM_E=1.0, erg_Guess=nothing, max_n_qnm=5)
    
    OmegaH = a ./ (2 .* (GNew .* M) .* (1 .+ sqrt.(1 .- a.^2)))
    alph = mu * GNew * M

    
    if (ergL(n, l, m, mu, M, a) < m .* OmegaH)||(for_s_rates==true)||QNM
        

        if (alph < 0.03)
            alph_ev = 0.03
        else
            alph_ev = alph
        end
        
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

        rescale = (alph ./ alph_ev).^(4 .* l + 5)
        
        b = sqrt.(1 - a.^2)
        
        
        if !QNM
#            if (n <= 5)&&(n >= 1)&&(l==0)&&(m==0)&&(alph > 0.03)&&(alph < 2.0)
#                erg0_r = open(readdlm, "input_info/ErgEvol_n_$(n)_l_0_m_0_.dat")
#                erg0_i = open(readdlm, "input_info/ErgEvol_I_n_$(n)_l_0_m_0_.dat");
#                itp = LinearInterpolation(erg0_r[:, 1], erg0_r[:, 2], extrapolation_bc=Line())
#                itpI = LinearInterpolation(erg0_i[:, 1], erg0_i[:, 2], extrapolation_bc=Line())
#                w0 =  itp(alph) .+ im .* itpI(alph)
#            else
            if isnothing(erg_Guess)
                SR211_g = sr_rates(n, l, m, alph_ev ./ (GNew * M), M, a)
                # print("test \t ", SR211_g, "\n")
                if SR211_g == 0
                    SR211_g = -1e-10 .* mu
                end
                w0 = (ergL(n, l, m, alph_ev ./ (GNew * M), M, a; full=false) .+ im * SR211_g) .* GNew * M
            else
                w0 = erg_Guess
            end
#            end
        else
            ### only valid for l= 0,m=0, otherwise better guesses are needed
            w0 = 0.117 .+ 0.004 .* (alph ./ 0.1) .- 0.084 * im .- 0.01 .* (alph ./ 0.1)  # rough guess for n = 1, l=0, m=0
        end
        
        
        
        function wrapper!(F, x)
            # wR = 10 .^ x[1]  # real
            wR = x[1]
            wI = x[2]  # imag
            erg = wR + im * wI
    
            q = - sqrt.(alph_ev.^2 .- erg.^2)
            
            if !QNM
                if real(q) > 0
                    q *= -1
                end
            else
                if real(q) < 0
                    q *= -1
                end
            end
           
            
            gam = im * a * sqrt.(erg.^2 .- alph_ev.^2)
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
            w0 = alph_ev + imag(w0)
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
        end
        
        if ((!sol.x_converged) && (!sol.f_converged))||(sol.zero[1] .< 0.0)
            if return_both
                return 0.0, 0.0
            else
                return 0.0
            end
        end

        if !return_both
            return sol.zero[2] .* rescale  # im(G * M * omega)
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
    
    
    alph = mu * GNew * M
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
        

        #### TESTING
        trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
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



function gf_radial(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=1000, Npts_Bnd=1000, debug=false, Ntot_safe=5000,  iter=10, xtol=1e-10, ftol=1e-10, tag="_", Nang=100000, eps_fac = 1e-10, m=0, l=0, NON_REL=false, h_mve=100, to_inf=false)
    
    rp = BigFloat(1.0 .+ sqrt.(1.0 .- a.^2))
    rmm = BigFloat(1.0 .- sqrt.(1.0 .- a.^2))
    alph = mu * GNew * M
    
    
    erg_1G = ergL(n1, l1, m1, mu, M, a; full=false)
    erg_2G = ergL(n2, l2, m2, mu, M, a; full=false)
    erg_3G = ergL(n3, l3, m3, mu, M, a; full=false)
    erg_pxy = erg_1G + erg_2G - erg_3G
    
    maxN = maximum([n1 n2 n3])
    minN = maximum([n1 n2 n3])
    
    
    ### quick check if to inf to not!
    erg1R, erg1I = find_im_part(mu, M, a, n1, l1, m1; Ntot_force=5000, for_s_rates=true, return_both=true)
    erg2R, erg2I = find_im_part(mu, M, a, n2, l2, m2; Ntot_force=5000, for_s_rates=true, return_both=true)
    erg3R, erg3I = find_im_part(mu, M, a, n3, l3, m3; Ntot_force=5000, for_s_rates=true, return_both=true)
    if (erg1R .+ erg2R .- erg3R) .> alph
        to_inf = true
    else
        to_inf = false
    end
    
    #### fix l,m
    if to_inf
        m = (m1 + m2 - m3)
        l = l1 + l2 - l3
        rmax = Float64.(100 ./ alph.^2 .* (minN ./ 2.0) ) .* 100.0
    else
        rmax = Float64.(100 ./ alph.^2 .* (minN ./ 2.0) ) 
    end
    
    
    rmax = Float64.(100 ./ alph.^2 .* (minN ./ 2.0) )
    
    rlist = 10 .^range(log10(rp .* (1.0 .+ eps_fac)), log10.(rmax), rpts)
    
    if NON_REL
        rf_1 = radial_bound_NR(n1, l1, m1, mu, M, rlist)
        erg_1 = erg_1G * GNew * M
    else
        rl, r1, erg_1 = solve_radial(mu, M, a, n1, l1, m1; rpts=Npts_Bnd, return_erg=true, Ntot_safe=Ntot_safe)
        itp = LinearInterpolation(log10.(rl), r1, extrapolation_bc=Line())
        rf_1 = itp(log10.(rlist))
        if imag(erg_1) < 0
            return 0.0
        end
    end
    
    if NON_REL
        rf_2 = radial_bound_NR(n2, l2, m2, mu, M, rlist)
        erg_2 = erg_2G * GNew * M
    else
        if (n2 == n1)&&(l2 == l1)&&(m2 == m1)
            rf_2 = rf_1
            erg_2 = erg_1
        else
            rl, r2, erg_2 = solve_radial(mu, M, a, n2, l2, m2; rpts=Npts_Bnd,  return_erg=true, Ntot_safe=Ntot_safe)
            itp = LinearInterpolation(log10.(rl), r2, extrapolation_bc=Line())
            rf_2 = itp(log10.(rlist))
            if imag(erg_2) < 0
                return 0.0
            end
        end
    end
    
    if NON_REL
        rf_3 = radial_bound_NR(n3, l3, m3, mu, M, rlist)
        erg_3 = erg_3G * GNew * M
    else
        rl, r3, erg_3 = solve_radial(mu, M, a, n3, l3, m3; rpts=Npts_Bnd, return_erg=true, Ntot_safe=Ntot_safe)
        itp = LinearInterpolation(log10.(rl), r3, extrapolation_bc=Line())
        rf_3 = itp(log10.(rlist))
        if imag(erg_3) < 0
            return 0.0
        end
    end
    
    erg = (erg_1 + erg_2 - erg_3) + 0 * im # leave the 0 im for NR case
    
    
    Z1 = spheroidals(l1, m1, a, erg_1)
    Z2 = spheroidals(l2, m2, a, erg_2)
    Z3 = spheroidals(l3, m3, a, erg_3)
    Z4 = spheroidals(l, m, a, erg)
   
    thetaV = acos.(1.0 .- 2.0 .* rand(Nang))
    phiV = rand(Nang) .* 2*pi

    function func_ang(x)
        return real(Z1.(x[1], x[2]) .* Z2.(x[1], x[2]) .* conj(Z3.(x[1], x[2])) .* conj(Z4.(x[1], x[2])))
    end
    function func_ang_2(x)
        return real(Z1.(x[1], x[2]) .* Z2.(x[1], x[2]) .* conj(Z3.(x[1], x[2])) .* conj(Z4.(x[1], x[2]))) .* a.^2 .* cos.(x[1]).^2
    end
    CG = 0.0
    CG_2 = 0.0
    for i in 1:Nang
        CG += func_ang([thetaV[i], phiV[i]])
        CG_2 += func_ang_2([thetaV[i], phiV[i]])
    end
    CG *= 4*pi / Nang
    CG_2 *= 4*pi / Nang
    
    
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
    # rout_star = r_list_map .+ 2.0 .* rp ./ (rp .- rmm) .* log.(r_list_map ./ rp .- 1.0) .- 2.0 .* rmm .* log.(r_list_map ./ rmm .- 1.0) ./ (rp .- rmm)
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
#        if add_source
#            newV += - net_rescale .* SS .* (CG .* r_input.^2 .+ CG_2 .* a.^2) .* delt ./ ff.^(3/2)
#        end
        
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
    
    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
    
    midP = Int(round(length(rvals) - 10))
    wronk  = (outWF_fw[midP] .* (outWF[midP+1] .- outWF[midP-1]) .-  outWF[midP] .* (outWF_fw[midP+1] .- outWF_fw[midP-1]) )./ (itp_rrstar.(rvals[midP+1]) .- itp_rrstar.(rvals[midP-1]))
    wronk *= (itp_rrstar.(rvals[midP]).^2 .- 2 .* itp_rrstar.(rvals[midP]) .+ a.^2) .* (GNew .* M)

    Tmm = (itpG(log10.(itp_rrstar.(rvals))) + im * itpGI(log10.(itp_rrstar.(rvals)))) .* (CG .* itp_rrstar.(rvals).^2 .+ CG_2 .* a.^2)
    
#    ###
#    wronk_test = []
#    r_test = []
#    for i in 2:(length(rvals)-1)
#        tt  = (outWF_fw[i] .* (outWF[i+1] .- outWF[i-1]) .-  outWF[i] .* (outWF_fw[i+1] .- outWF_fw[i-1]) )./ (2 * (itp_rrstar.(rvals[i+1]) .- itp_rrstar.(rvals[i-1])) )
#        tt *= (itp_rrstar.(rvals[i]).^2 .- 2 .* itp_rrstar.(rvals[i]) .+ a.^2) .* (GNew .* M)
#        append!(wronk_test, tt)
#        append!(r_test, itp_rrstar.(rvals[i]))
#    end
#    writedlm("test_store/testW_r.dat",hcat(r_test, real(wronk_test)))
#    writedlm("test_store/testW_i.dat",hcat(r_test, real(wronk_test)))
#    ##
    
    ####
    if to_inf
        
        out_R = outWF[end] .* trapz(outWF_fw .* Tmm, itp_rrstar.(rvals)) ./ wronk
        maxV = Float64.(real(out_R .* conj.(out_R)))

        lam = (mu ./ (M_pl .* 1e9))^2
        kk = real(sqrt.(erg.^2 .- alph.^2))
        rate_out = 2 .* alph .* kk .* (maxV[end] .* itp_rrstar.(itp_rrstar.(rvals[end])).^2) .* lam^2
        out_gamma = rate_out ./ mu^2 .* (GNew * M^2 * M_to_eV)^2
        
        ### save full WF?
#        outR = []
#        r_out = []
#        nskip = Int(round(length(rvals) / 100))
#        println("filling output ")
#        for i in 2:nskip:(length(rvals)-1)
#            temp = (outWF_fw[i] .* trapz(outWF[i:end] .* Tmm[i:end], itp_rrstar.(rvals[i:end])) .+ outWF[i] .* trapz(outWF_fw[1:i] .* Tmm[1:i], itp_rrstar.(rvals[1:i]))) ./ wronk
#            append!(outR, temp)
#            append!(r_out, rvals[i])
#        end
#        println("done filling output ")
#        maxV = Float64.(real(outR .* conj.(outR)))
#        writedlm("test_store/test_1.dat", hcat(itp_rrstar.(r_out), maxV))
        ####
    else
        # idx_hold = itp_rrstar.(rvals) .> 1.01 .* rp
        # rnew_rp = trapz(outWF[idx_hold] .* Tmm[idx_hold] , itp_rrstar.(rvals[idx_hold])) .* outWF_fw[idx_hold][1] ./ wronk
        rnew_rp = trapz(outWF .* Tmm , itp_rrstar.(rvals)) .* outWF_fw[1] ./ wronk

        maxV = real(rnew_rp.* conj.(rnew_rp))

        lam = (mu ./ (M_pl .* 1e9))^2
        rate_out = 4 .* alph.^2 .* (1 .+ sqrt.(1 - a.^2)) .* Float64(maxV) .* lam^2
        out_gamma = rate_out ./ mu^2 .* (GNew * M^2 * M_to_eV)^2
        if debug
            # theory-NR
            if (n1==2)&&(n2==2)&&(n3==3) && (l1==1)&&(l2==1)&&(l3==2)
                test = 4.3e-7 * alph.^(11) .* (1 .+ sqrt.(1 .- a.^2)) # good
            elseif (n1==2)&&(n2==4)&&(n3==3) && (l1==1)&&(l2==1)&&(l3==2)
                test = 2.5e-8 * alph.^(11) .* (1 .+ sqrt.(1 .- a.^2)) # Don't agree... smaller? by O(10)? NR lim, inf dominated!
            elseif (n1==4)&&(n2==4)&&(n3==3) && (l1==1)&&(l2==1)&&(l3==2)
                test = 1.7e-11 * alph.^(11) .* (1 .+ sqrt.(1 .- a.^2)) ##  Not agreeing, but i dont agree with Masha? NR lim, inf comparable?
            elseif (n1==2)&&(n2==2)&&(n3==4) && (l1==1)&&(l2==1)&&(l3==2)
                test = 1.5e-7 * alph.^(11) .* (1 .+ sqrt.(1 .- a.^2)) ## significantly smaller? Bound dominated
            elseif (n1==4)&&(n2==4)&&(n3==4) && (l1==1)&&(l2==1)&&(l3==2)
                test = 2.3e-7 * alph.^(7) .* (1 .+ sqrt.(1 .- a.^2)) ## significantly smaller? Resonance...
            elseif (n1==2)&&(n2==3)&&(n3==4) && (l1==1)&&(l2==2)&&(l3==3)
                test = 9.1e-8 * alph.^(11) .* (1 .+ sqrt.(1 .- a.^2)) ## Bound dominated
            elseif (n1==2)&&(n2==4)&&(n3==4) && (l1==1)&&(l2==2)&&(l3==3)
                test = 7.8e-11 * alph.^7 .* (1 .+ sqrt.(1 .- a.^2)) # Resonance..., small by O(10)
            else
                test = 0.0
            end
            print("NR Rate: ", test, "\n")
        end
    end

    if debug
        print("Rate \t", out_gamma, "\n")
    end
    return out_gamma # unitless [gamma / mu]
end

