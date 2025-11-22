include("solve_sr_rates.jl")

rr = [(2,1,1), (3,1,1), (3,2,2)]
debug = true
xtol=1e-15
ftol=1e-30
iter=100

for i in 1:length(rr)
    alphaL = LinRange(0.05, 0.2, 10)
    M = 137;
    M2 = 100;
    a = 0.9;
    rp = 1.0 .+ sqrt.(1 - a .^2)
    n, l, m = rr[i]
    Ntot_safe = 6000
    rpts=7000

    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
    tempH = zeros(length(alphaL), 3)
    
    for j in 1:length(alphaL)
    	rl, r1, erg = solve_radial(alphaL[j] ./ (GNew .* M), M, a, n, l, m; rpts=rpts, return_erg=true, Ntot_safe=Ntot_safe)
        wR = real.(erg)
        wI = imag.(erg)
        rl2, r12, erg2 = solve_radial(alphaL[j] ./ (GNew .* M), M, a, n, l, m-2; rpts=rpts, return_erg=true, Ntot_safe=Ntot_safe)
        wR2 = real.(erg2)
        wI2 = imag.(erg2)
        
        period_res = 2 * pi ./ (abs.(wR .- wR2) ./ (GNew .* M))
        r_res = ((period_res ./ (2 * pi)).^2 .* GNew .* (M .+ M2) ) .^(1/3) ./ (GNew .* M)

        factor1 = sqrt.(abs.(wI2 ./ wI))
        ergD = wR2 .- wR
        println("here \t", wR)
        
        itp_base = interpolate((log10.(rl),), r1, Gridded(Linear()))
        itp = extrapolate(itp_base, 0.0)  # Zero outside the domain
        
        itp_base2 = interpolate((log10.(rl2),), r12, Gridded(Linear()))
        itp2 = extrapolate(itp_base2, 0.0)  # Zero outside the domain
        
        
        ## SOLVE....
        function wrapper!(F, x)
            # println(x)
            if (i == 1)||(i==2)
                af = 0.31
            elseif i == 3
                af = -0.18
            end
            
            if x[1] < log10.(2 .* rp)
                x[1] = log10.(2 * rp)
            end
            rlist = 10 .^range(log10.(rp .* (1.0 .+ 1.1)), log10.(x), 2000)
            
            rf_1 = itp(log10.(rlist))
            rf_1[rlist .> maximum(rl)] .= 0.0
            rf_2 = itp2(log10.(rlist))
            rf_2[rlist .> maximum(rl2)] .= 0.0
            rf_1[isnan.(rf_1)] .= 0.0
            rf_2[isnan.(rf_2)] .= 0.0
            
            # println(rf_1)
            # println(rf_2)
            
            factor2 = trapz(rf_1 .* conj.(rf_2) .* rlist.^2 .* (rlist.^2 ./ x.^(2+1)), rlist)
            potentialV = af .* factor2 .* sqrt.(6 * pi / 5) / 2 .* alphaL[j]
            delwI = abs.((wI .- wI2) ./ (wR .- wR2).^2 .* abs.(potentialV).^2)
            if !isnan.(potentialV)
                F[1] = abs.(wI .- delwI) ./ abs.(wI)
            else
                F[1] = 1
            end
            # println("x/f \t",x, "\t", F, "\t", factor2)
        end
        
        # r_test = 10 .^range(log10.(2 .* rp .* (1.0 .+ 1.1)), log10.(2e5), 100)
        r_test = 10 .^range(log10.(1e3), log10.(2e5), 10000)
        out_test = zeros(length(r_test))
        for i in 1:length(r_test)
            F = [0.0]
            wrapper!(F, r_test[i])
            out_test[i] = F[1]
        end
        full_out = [r_test[argmin(out_test)]]
        
        
        # println(maximum(rl))
        # sol = nlsolve(wrapper!, [log10.(maximum(rl) ./ 1.2) ], show_trace=false, autodiff = :forward, xtol=xtol, ftol=ftol, iterations=iter)
        # full_out = 10 .^ sol.zero
        println([alphaL[j], full_out[1], r_res])
        tempH[j,:] .= [alphaL[j], full_out[1], r_res] # everything in (G*M) units
    end
    writedlm("test_store/binary_test_$(i).dat", tempH)
end
