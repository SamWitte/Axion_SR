include("solve_sr_rates.jl")
using DelimitedFiles
using SpecialFunctions
using HypergeometricFunctions
using Interpolations
using PyCall
sp = pyimport("scipy.special")

mu = 2.5e-13
M = 22.2
a = 0.9

n1 = 2
l1 = 1
m1 = 1

n2 = 2
l2 = 1
m2 = 1

n3 = 3
l3 = 2
m3 = 2

kpts=40
rpts=2000
rmaxT=100
include_cont = false
nmax=4

function dRbnd(n, mu, M, rlist; first=true)
    # assume l = 0 m = 0
    alph = M * GNew * mu
    a0 = 1 ./ alph.^2
    if first
        h1 = sp.hyp1f1(1 - n, 2, 2 .* rlist ./ (n .* a0))
        h2 = sp.hyp1f1(2 - n, 3, 2 .* rlist ./ (n .* a0))
        out = exp.(-rlist ./ (n .* a0)) .* (-2 .* h1 .- 2 .* (n .- 1) .* h1) ./ (n .* a0).^(5/2)
        return out
    else
        _genlaguerre(n, α, x) = binomial(n+α,n) * HypergeometricFunctions.M.(-n, α+1, x)
        h1 = sp.hyp1f1(2 - n, 4, 2 .* rlist ./ (n .* a0))
        h2 = _genlaguerre(-1 + n, 1, 2 .* rlist ./ (n .* a0))
        out = 2 .* exp.(-rlist ./ (n .* a0)) .* (2/3 .* n .* (n.^2 .- 1) .* h1 .+ h2) ./ (n.^(9/2) .* a0.^(7/2))
        return out
    end
end

function dRcon(rlist, baseF; first=true)
    
    itp = LinearInterpolation(log10.(rlist), baseF, extrapolation_bc=Line())
    rnew = LinRange(rlist[1], rlist[end], 1000000)
    outOld = itp(log10.(rnew))
    hh = rnew[2] .- rnew[1]
    outNew = zeros(Complex, length(rnew))
    if first
        for i in 2:(length(rnew)-1)
            outNew[i] = (outOld[i+1] .- outOld[i-1]) ./ (2 * hh)
        end
        outNew[1] = outNew[2]
        outNew[end] = outNew[end-1]
    else
        for i in 2:(length(rnew)-1)
            outNew[i] = (outOld[i+1] .- 2 .* outOld[i] .+ outOld[i-1]) ./ hh.^2
        end
        outNew[1] = outNew[2]
        outNew[end] = outNew[end-1]
    end
    
    itp = LinearInterpolation(rnew, outNew, extrapolation_bc=Line())
    return itp(rlist)
end

function test_radial_reconst(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; kpts=10, rpts=2000, rmaxT=100, inf_nr=true, Nang=100000, Npts_Bnd=1000, debug=false, include_cont=true, Ntot_safe=5000, sve_for_test=false, bnd_thresh=1e-3, use_analytic=false, eps_r=1e-10, nmax=10)
    
    rp = 1 + sqrt.(1 - a.^2)
    alph = mu * GNew * M
    print("alpha \t", alph, "\n")
    
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
    
    rf_1 = radial_bound_NR(n1, l1, m1, mu, M, rlist)
    erg_1 = erg_1G .* GNew * M
    
    rf_2 = radial_bound_NR(n2, l2, m2, mu, M, rlist)
    erg_2 = erg_2G .* GNew * M
    
    rf_3 = radial_bound_NR(n3, l3, m3, mu, M, rlist)
    erg_3 = erg_3G .* GNew * M
    
    
    erg_ind = erg_1 .+ erg_2 - erg_3
    k_ind_2 = alph.^2 .- erg_ind.^2
    
   

    Z1 = spheroidals(l1, m1, a, erg_1 ./ (GNew .* M))
    Z2 = spheroidals(l2, m2, a, erg_2 ./ (GNew .* M))
    Z3 = spheroidals(l3, m3, a, erg_3 ./ (GNew .* M))
    
    
    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
        

    # compute bound contribution
    kmin = 0.02 .* alph^2
    kmax = 2.0 .* alph^2
    # kk_list = 10 .^LinRange(log10.(kmin), log10.(kmax),  kpts)
    kk_list = LinRange(kmin, kmax,  kpts)
    dk = kk_list[2] .- kk_list[1]
    
    ck_list = zeros(Complex, kpts)
 
    UnF = GNew .* M
    bound_c = 0.0
    fullWF = zeros(length(rlist))
    
    LHS_check =  zeros(length(rlist))
    contin_c = 0.0
    done_nmax = false

    
    function func_ang(x, Zf1, Zf2, Zf3, Zf4)
        return real(Zf1.(x[1], x[2]) .* Zf2.(x[1], x[2]) .* conj(Zf3.(x[1], x[2])) .* conj(Zf4.(x[1], x[2])))
    end
        
        
    for n in 1:nmax
        
        rf_4 = radial_bound_NR(n, 0, 0, mu, M, rlist)
        
        
        erg_4G = ergL(n, 0, 0, mu, M, a)
        erg_4 = erg_4G .* GNew * M
        Z4 = spheroidals(0, 0, a, erg_4 ./ (GNew .* M))
        
        thetaV = acos.(1.0 .- 2.0 .* rand(Nang))
        phiV = rand(Nang) .* 2*pi
        CG = 0.0
        for i in 1:Nang
            CG += func_ang([thetaV[i], phiV[i]], Z1, Z2, Z3, Z4)
        end
        CG *= 4*pi / Nang
       
        
        r_integrd = (rf_1 ./ UnF^(3/2)) .* (rf_2 / UnF^(3/2)) .* conj(rf_3 / UnF^(3/2)) .* conj(rf_4 / UnF^(3/2)) .* (rlist * UnF).^2
        radial_int = trapz(r_integrd, rlist * UnF)
        
        
        kdiff_sq = (erg_4.^2 - erg_ind^2) ./ (GNew .* M).^2
        
        ff = CG .* radial_int ./ (2 .* mu).^(3 / 2) ./ (kdiff_sq)
        
        SS = CG ./ (2 .* alph).^(3/2) .* rf_1 .* rf_2 .* conj(rf_3)
        
        ### double factor 1/2
        if (n1==n2)&&(l1==l2)&&(m1==m2)
            ff /= 2
            SS /= 2
        end
        
        
        res_n = (rf_4[1] ./ UnF^(3/2)) .* ff
        if debug
            print(n, "\t", Float64(abs.(res_n)), "\n")
        end
        
        bound_c += res_n
        
        
        
        
        kdiff_sq2 = (erg_4.^2 - erg_ind.^2)
        r_integrd2 = rf_1 .* rf_2 .* conj(rf_3) .* conj(rf_4) .* rlist.^2
        radial_int2 = trapz(r_integrd2, rlist)
        ff2 = CG .* radial_int2 ./ (2 .* alph).^(3 / 2) ./ (kdiff_sq2)
        Rad = rf_4 .* ff2
        
        fullWF += Rad
        
        LHS_check .+= (alph.^2 .- erg_ind.^2) .* Rad
        LHS_check .+= - 2 .* alph.^2 .* Rad ./ rlist
        LHS_check .+= - 2 .* dRbnd(n, mu, M, rlist; first=true) ./ rlist .* ff2
        LHS_check .+= - dRbnd(n, mu, M, rlist; first=false) .* ff2
        
        
        
        if abs.(res_n / bound_c) < bnd_thresh
            done_nmax = true
        end
        
        if n == 1
            writedlm("test_store/Srs_test.dat", hcat(rlist, SS))
        end
        print("N lvl \t", n, "\n")
        
    end
    


    if include_cont
        # compute continuous contribution
        for i in 1:kpts
            k = kk_list[i]
            erg_New = sqrt.(k.^2 .+ alph.^2)
            Z4 = spheroidals(0, 0, a, erg_New ./ (GNew .* M))
            


            thetaV = acos.(1.0 .- 2.0 .* rand(Nang))
            phiV = rand(Nang) .* 2*pi

            CG = 0.0
            
            for i in 1:Nang
                CG += func_ang([thetaV[i], phiV[i]], Z1, Z2, Z3, Z4)
            end
            CG *= 4*pi / Nang
        
            
            out_goingR = radial_inf_NR(k, 0, mu, M, rlist)
            
            r_integrd = rf_1 .* rf_2 .* conj(rf_3) .* conj(out_goingR) .* rlist.^2
            radial_int = trapz(r_integrd, rlist)
            
            
            ff = CG .* radial_int ./ (2 .* alph).^(3 / 2) ./ (k_ind_2 .+ k.^2) ./ (2 * pi)
            
            ### double factor 1/2
            if (n1==n2)&&(l1==l2)&&(m1==m2)
                ff /= 2
            end
            
            if !isnan.(out_goingR[1] .* ff)
                ck_list[i] = out_goingR[1] .* ff
            else
                print("getting nan... \t", i, "\t", out_goingR[1:3], "\t", ff, "\n")
            end
            if debug
                print("cont \t", i, "\t", kk_list[i], "\t", Float64(abs.(ck_list[i])), "\t", ff, "\n")
            end
            
            Rad = dk .* ff .* out_goingR
            LHS_check .+= real((alph.^2 .- erg_ind.^2) .* Rad)
            LHS_check .+= - real(2 .* alph.^2 .* Rad ./ rlist)
            
            
            LHS_check .+= - real(2 .* dRcon(rlist, Rad; first=true) ./ rlist)
            LHS_check .+= - real(dRcon(rlist, Rad; first=false) )
            
        end
        contin_c += trapz(ck_list, kk_list)
    end
    
    psi_1 = bound_c .+ contin_c
    print("bound contribution:  ", Float64(abs.(bound_c)), "\t contin: ", Float64(abs.(contin_c)), "\n")
        
    
    lam = (mu ./ (M_pl .* 1e9))^2
    rate_out = 4 .* alph.^2 .* (1 .+ sqrt.(1 - a.^2)) .* Float64(real(psi_1 .* conj(psi_1))) .* lam^2
    # rate_out = 2 * alph.^2 .* rp * Float64(real(psi_1 .* conj(psi_1))) .* lam^2
    print("alph / rate \t", alph, "\t", rate_out ./ mu^2 .* (GNew * M^2 * M_to_eV)^2, "\n")
    print("mu2-erg2 \t", alph.^2 .- erg_ind.^2, "\n")
    
    ### test if it works....
    
    return fullWF ./ (GNew .* M), LHS_check, rlist # unitless [gamma / mu]
end


wfO, LHS_check, rlist = test_radial_reconst(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; kpts=kpts, rpts=rpts, rmaxT=rmaxT, include_cont=include_cont, nmax=nmax)

writedlm("test_store/AttemptWF.dat", hcat(rlist, wfO))
writedlm("test_store/AttemptWF_err.dat", hcat(rlist, LHS_check))
