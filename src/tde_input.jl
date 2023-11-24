using Printf
using Roots
using SpecialFunctions
using DelimitedFiles
using PyCall
using Plots
using Dates
numpy = pyimport("numpy")

const Ms = 1.9891e30
const Rs = 6.957e8
const c = 2.998e8
const G = 6.674e-11



function tde_like(MassBH, aBH; plot=false)
    N_psi = 100  # Number of angles I throw stars in at (linearly spaced between 0 and pi/2).
    N_spin = 150  # Number of absolute values of spins linearly spaced between 0 and 0.9999
    
    
    nameF = "input_info/max_mass_matrix_$(N_psi)_$(N_spin).txt"
    if isfile(nameF)
        max_mass_matrix = readdlm(nameF)
    else
        _, max_mass_matrix = get_hills_masses(N_psi, N_spin)  # Need this matrix for tidal force likelihood
        writedlm(nameF, max_mass_matrix)
    end

    prior_spins(a) = ones(size(a))  # agnostic spin prior.
    if plot
        p_a, a = one_d_spin_fixed_mass(1e8 * Ms, prior_spins, max_mass_matrix)
        p = plot(a[p_a .> 0], p_a[p_a .> 0], label="1e8")
        p_a, a = one_d_spin_fixed_mass(2e8 * Ms, prior_spins, max_mass_matrix)
        plot!(a[p_a .> 0], p_a[p_a .> 0], label="2e8")
        p_a, a = one_d_spin_fixed_mass(3e8 * Ms, prior_spins, max_mass_matrix)
        plot!(a[p_a .> 0], p_a[p_a .> 0], label="3e8")
        p_a, a = one_d_spin_fixed_mass(4e8 * Ms, prior_spins, max_mass_matrix)
        plot!(a[p_a .> 0], p_a[p_a .> 0], label="4e8")
        xlims!((0, 1))
        ylims!((0, 20))
        savefig("/Users/samuelwitte/Desktop/Test.png")
    end
    
    p_a, a = one_d_spin_fixed_mass(MassBH * Ms, prior_spins, max_mass_matrix)

    prob = p_a[argmin(abs.(a .- aBH))]
    if prob <= 0
        prob = 1e-100
    end
    return log.(prob)
end

function get_ibso(a, psi)
    cs = [a^8 * cos(psi)^4,
          0,
          -2 * a^6 * cos(psi)^2 + 6 * a^6 * cos(psi)^4,
          -8 * a^4 * cos(psi)^2 + 16 * a^4 * cos(psi)^2 * sin(psi)^2,
          a^4 - 4 * a^4 * cos(psi)^2 + 9 * a^4 * cos(psi)^4,
          8 * a^2 - 24 * a^2 * cos(psi)^2 - 16 * a^2 * sin(psi)^2,
          16 - 2 * a^2 + 6 * a^2 * cos(psi)^2,
          -8,
          1]

    root = numpy.roots(reverse(cs))
    

    if a > 1
        return NaN
    end

    real_roots = real.(root[isreal.(root)])
   
    if length(real_roots) != 2
        if a > 0
            return real_roots[2]
        end
        return real_roots[1]
    end

    if a > 0
        return minimum(real_roots)
    end
    return maximum(real_roots)
end

function hills_mass(a, psi, eta=1, mstar=1, rstar=1)
    x = get_ibso(a, psi)
    amp = (2/eta * c^6 * Rs^3 / (G^3 * Ms)) ^ 0.5 * 1 / x .^1.5
    fac = (1 + 6 * x ./ (x .^2 - a^2 * cos(psi)^2) + 3 * a^2 / (2 * x .^2) - 6 * a * sin(psi) / x .^1.5 * (x .^2 / (x .^2 - a^2 * cos(psi)^2)) .^ 0.5) .^ 0.5
    
    return amp * fac
end

function get_hills_masses(N_spin=300, N_psi=200, Mstar=1)
    psi_ = range(0.001, stop=pi/2, length=N_psi)
    a_ = range(-0.9999, stop=0.9999, length=2 * N_spin)
    Mstar *= Ms

    @printf("Generating Hills masses......\n")
    max_mass_matrix = zeros(length(a_), length(psi_))

    for i in 1:length(a_)
        for j in 1:length(psi_)
            max_mass_matrix[i, j] = hills_mass(a_[i], psi_[j])
        end
    end

    max_m = zeros(size(a_))
    for j in 1:length(a_)
        max_m[j] = maximum(max_mass_matrix[j, :])
    end

    f = (MassRadiusRelation(Mstar) ./ Rs) .^ 1.5 * (Ms ./ Mstar) .^ 0.5
    max_m .*= f

    return max_m, max_mass_matrix
end

function MassRadiusRelation(mstar)
    rstar = zeros(size(mstar))
    
    if sum(mstar .< Ms) > 0
        rstar[mstar .< Ms] .= Rs * (mstar[mstar .< Ms] / Ms) .^ 0.56
    end
    if sum(mstar .> Ms) > 0
        rstar[mstar .> Ms] .= Rs * (mstar[mstar .> Ms] / Ms) .^ 0.79
    end
    return rstar
end

function KroupaIMF(mstar)
    i_lowest = mstar .< 0.08 * Ms
    i_low = (mstar .< 0.5 * Ms) .& .! i_lowest
    i_med = (mstar .< 1 * Ms) .& .! i_lowest .& .! i_low
    i_high = .! i_lowest .& .! i_low .& .! i_med

    P = zeros(size(mstar))

    P[i_lowest] .= (mstar[i_lowest] / (0.08 * Ms)) .^ (-0.3)
    P[i_low] .= (mstar[i_low] / (0.08 * Ms)) .^ (-1.8)
    P[i_med] .= (mstar[i_med] / (0.08 * Ms)) .^ (-2.7) * (0.5 / 0.08) ^ (-1.8 + 2.7)
    P[i_high] .= (mstar[i_high] / (0.08 * Ms)) .^ (-2.3) * (1 / 0.08) ^ (-2.7 + 2.3) * (0.5 / 0.08) ^ (-1.8 + 2.7)

    return P
end

function TDE_rate_on_stellar_params(mstar)
    rstar = MassRadiusRelation(mstar)
    return (mstar / Ms) .^ (-1 / 3) .* (rstar / Rs) .^ 0.25
end

function StellarMassDistributionFunction(;Mmin=0.01 * Ms, Mmax=2.0 * Ms, N=10000)
    mstars = 10 .^LinRange(log10(Mmin), log10(Mmax), N)
    dM = mstars[2:end] - mstars[1:end-1]
    ps = KroupaIMF(mstars) .* TDE_rate_on_stellar_params(mstars)
    return ps[1:end-1] ./ sum(ps[1:end-1] .* dM)
end

function one_d_spin_fixed_mass(Mbh, prior_spin, max_mass_matrix; prior_star=nothing, prior_psi=nothing,
                                N_star=100, Mstar_min=0.1, Mstar_max=1, eta=1)

    Mstar_min *= Ms
    Mstar_max *= Ms

    mstars = 10 .^LinRange(log10(Mstar_min), log10(Mstar_max), N_star + 1)
    
    if isnothing(prior_star)
        psd = StellarMassDistributionFunction(Mmin=Mstar_min, Mmax=Mstar_max, N=N_star + 1)
    else
        psd = prior_star(mstars)
    end

    dMs = mstars[2:end] - mstars[1:end-1]
    mstars = mstars[1:end-1]

    N_spin = size(max_mass_matrix, 1)
    N_psi = size(max_mass_matrix, 2)

    psi_ = LinRange(0.001, pi/2, N_psi)
    a_ = LinRange(-0.9999, 0.9999, N_spin)

    dpsi = psi_[2] - psi_[1]
    da_ = a_[2] - a_[1]

    f = (MassRadiusRelation(mstars) ./ Rs) .^ 1.5 .* (Ms ./ mstars) .^ 0.5 .* eta .^ 0.5
    
    
    
    
    mass_mass_matrix = []
    for k in 1:N_spin
        tmp_matrix = zeros(N_psi, N_star)
        for u in 1:N_psi
            tmp_matrix[u, :] .= max_mass_matrix[k, u] * f
        end
        push!(mass_mass_matrix, tmp_matrix)
    end

    p_a = zeros(N_spin)

    time1=Dates.now()
    
    for i in 1:N_spin
        
        p_ = zeros(size(mass_mass_matrix[i]))
        p_[mass_mass_matrix[i] .- Mbh .> 0.0] .= 1.0
        val = sum(sum(p_ .* transpose(psd[:]) .* transpose(dMs[:]), dims=2)[:] .* cos.(psi_) .* dpsi, dims=1) .* prior_spin(a_[i]) .* da_
        
        p_a[i] = val[1]
    end
    # print(size(p_a[Int(N_spin / 2)+1:end]),"\t", size(reverse(p_a[1:Int(N_spin / 2)])), "\n")
    p_a = p_a[Int(N_spin / 2)+1:end] .+ reverse(p_a[1:Int(N_spin / 2)])
    
    
    p_a ./= sum(p_a .* da_)

    return p_a, a_[Int(N_spin / 2)+1:end]
end



# tde_like(1e8, .5)
