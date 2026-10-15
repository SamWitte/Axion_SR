# Pointwise linear operators acting on a metric perturbation h_{mu nu}(r, th) e^{-i w v + i m phi}, written as
# coefficients on the 60 "jets" (10 components x {value, d_r, d_th, d_rr, d_rth, d_thth}).
#   Component order c: vv, vr, vth, vph, rr, rth, rph, thth, thph, phph.

const PAIRS = [(1, 1), (1, 2), (1, 3), (1, 4), (2, 2), (2, 3), (2, 4), (3, 3), (3, 4), (4, 4)]
const CIDX = let M = zeros(Int, 4, 4)
    for (c, (i, j)) in enumerate(PAIRS)
        M[i, j] = c; M[j, i] = c
    end
    M
end
const NJ = 60
jet(c, d) = (c - 1) * 6 + d          # d: 1 val, 2 r, 3 th, 4 rr, 5 rth, 6 thth

"""Jet arrays: H[μ,ν,:], dH[ρ,μ,ν,:] = ∂_ρ h_μν, ddH[σ,ρ,μ,ν,:] = ∂_σ∂_ρ h_μν, for frequency w and azimuthal m."""
function jet_arrays(w, m)
    H = zeros(ComplexF64, 4, 4, NJ)
    dH = zeros(ComplexF64, 4, 4, 4, NJ)
    ddH = zeros(ComplexF64, 4, 4, 4, 4, NJ)
    kv = -im * w; kp = im * m          # ∂_v, ∂_phi
    for μ in 1:4, ν in 1:4
        c = CIDX[μ, ν]
        e(d) = (x = zeros(ComplexF64, NJ); x[jet(c, d)] = 1; x)
        H[μ, ν, :] = e(1)
        first = [kv * e(1), e(2), e(3), kp * e(1)]
        for ρ in 1:4
            dH[ρ, μ, ν, :] = first[ρ]
        end
        # second derivatives: indices (σ,ρ)
        sec = Dict((1, 1) => kv^2 * e(1), (1, 2) => kv * e(2), (1, 3) => kv * e(3), (1, 4) => kv * kp * e(1),
                   (2, 2) => e(4), (2, 3) => e(5), (2, 4) => kp * e(2),
                   (3, 3) => e(6), (3, 4) => kp * e(3), (4, 4) => kp^2 * e(1))
        for σ in 1:4, ρ in 1:4
            ddH[σ, ρ, μ, ν, :] = sec[(min(σ, ρ), max(σ, ρ))]
        end
    end
    return H, dH, ddH
end

"""First and second covariant derivatives of h (as jets): D1[ρ,μ,ν,:] = ∇_ρ h_μν, D2[σ,ρ,μ,ν,:] = ∇_σ∇_ρ h_μν."""
function cov_derivs(bg::BGPoint, H, dH, ddH)
    Γ = bg.Gam; dΓ = bg.dGam
    D1 = copy(dH)
    for ρ in 1:4, μ in 1:4, ν in 1:4, λ in 1:4
        c1 = Γ[λ, ρ, μ]; c2 = Γ[λ, ρ, ν]
        c1 != 0 && (D1[ρ, μ, ν, :] .-= c1 .* H[λ, ν, :])
        c2 != 0 && (D1[ρ, μ, ν, :] .-= c2 .* H[μ, λ, :])
    end
    D2 = copy(ddH)
    for σ in 1:4, ρ in 1:4, μ in 1:4, ν in 1:4
        acc = view(D2, σ, ρ, μ, ν, :)
        for λ in 1:4
            # ∂_σ of the Christoffel terms in ∇_ρ h_μν
            d1 = dΓ[σ][λ, ρ, μ]; d2 = dΓ[σ][λ, ρ, ν]
            d1 != 0 && (acc .-= d1 .* view(H, λ, ν, :))
            d2 != 0 && (acc .-= d2 .* view(H, μ, λ, :))
            g1 = Γ[λ, ρ, μ]; g2 = Γ[λ, ρ, ν]
            g1 != 0 && (acc .-= g1 .* view(dH, σ, λ, ν, :))
            g2 != 0 && (acc .-= g2 .* view(dH, σ, μ, λ, :))
            # connection acting on the 3-tensor ∇_ρ h_μν
            t1 = Γ[λ, σ, ρ]; t2 = Γ[λ, σ, μ]; t3 = Γ[λ, σ, ν]
            t1 != 0 && (acc .-= t1 .* view(D1, λ, μ, ν, :))
            t2 != 0 && (acc .-= t2 .* view(D1, ρ, λ, ν, :))
            t3 != 0 && (acc .-= t3 .* view(D1, ρ, μ, λ, :))
        end
    end
    return D1, D2
end

"""
Linearised Einstein tensor dG_μν[h] and the Lorenz-gauge-fixed operator
    E_μν = dG_μν - ∇_(μ Z_ν) + 1/2 g_μν ∇·Z,   Z_ν = ∇^μ hbar_μν      (= -1/2 □hbar_μν - R_μανβ hbar^αβ)
as jets.  Vacuum background (R_μν = 0).  Also returns Z_ν as jets.
"""
function einstein_jets(bg::BGPoint, D1, D2; kap=0.0, nlow=nothing)
    g = bg.g; gi = bg.gi
    dR = zeros(ComplexF64, 4, 4, NJ)
    for μ in 1:4, ν in 1:4, λ in 1:4, κ in 1:4
        gik = gi[λ, κ]
        gik == 0 && continue
        dR[μ, ν, :] .+= 0.5 .* gik .* (D2[λ, μ, κ, ν, :] .+ D2[λ, ν, κ, μ, :] .- D2[λ, κ, μ, ν, :] .- D2[ν, μ, κ, λ, :])
    end
    trR = zeros(ComplexF64, NJ)
    for μ in 1:4, ν in 1:4
        gi[μ, ν] != 0 && (trR .+= gi[μ, ν] .* dR[μ, ν, :])
    end
    dG = similar(dR)
    for μ in 1:4, ν in 1:4
        dG[μ, ν, :] = dR[μ, ν, :] .- 0.5 .* g[μ, ν] .* trR
    end
    # ∇_μ∇_κ hbar_αν = D2[μ,κ,α,ν] - 1/2 g_αν g^ρσ D2[μ,κ,ρ,σ]
    trD2 = zeros(ComplexF64, 4, 4, NJ)       # g^ρσ D2[μ,κ,ρ,σ]
    trD1 = zeros(ComplexF64, 4, NJ)          # g^ρσ D1[κ,ρ,σ] = ∇_κ h
    for μ in 1:4, κ in 1:4, ρ in 1:4, σ in 1:4
        gi[ρ, σ] != 0 && (trD2[μ, κ, :] .+= gi[ρ, σ] .* D2[μ, κ, ρ, σ, :])
    end
    for κ in 1:4, ρ in 1:4, σ in 1:4
        gi[ρ, σ] != 0 && (trD1[κ, :] .+= gi[ρ, σ] .* D1[κ, ρ, σ, :])
    end
    Z = zeros(ComplexF64, 4, NJ)             # Z_ν = g^{ακ} (∇_κ h_αν - 1/2 g_αν ∇_κ h)
    dZ = zeros(ComplexF64, 4, 4, NJ)         # ∇_μ Z_ν
    for ν in 1:4, α in 1:4, κ in 1:4
        gak = gi[α, κ]
        gak == 0 && continue
        Z[ν, :] .+= gak .* (D1[κ, α, ν, :] .- 0.5 .* g[α, ν] .* trD1[κ, :])
        for μ in 1:4
            dZ[μ, ν, :] .+= gak .* (D2[μ, κ, α, ν, :] .- 0.5 .* g[α, ν] .* trD2[μ, κ, :])
        end
    end
    divZ = zeros(ComplexF64, NJ)
    for μ in 1:4, ν in 1:4
        gi[μ, ν] != 0 && (divZ .+= gi[μ, ν] .* dZ[μ, ν, :])
    end
    E = similar(dR)
    for μ in 1:4, ν in 1:4
        E[μ, ν, :] = dG[μ, ν, :] .- 0.5 .* (dZ[μ, ν, :] .+ dZ[ν, μ, :]) .+ 0.5 .* g[μ, ν] .* divZ
    end
    # constraint damping (Gundlach et al. 2005, rho = -1): + kap n_(μ Z_ν), n future-directed unit normal (covector).
    # Vanishes on solutions with Z = 0; turns □Z = 0 into a damped wave equation, removing near-zero Z modes.
    if kap != 0
        for μ in 1:4, ν in 1:4
            E[μ, ν, :] .+= (kap / 2) .* (nlow[μ] .* Z[ν, :] .+ nlow[ν] .* Z[μ, :])
        end
    end
    return E, dG, Z
end

"""Coefficients (on the jets of h) of S = -δ□[h]φ = h^μν ∇_μ∇_νφ + Z^ν ∂_νφ, given dphi[μ] = ∂_μφ, ddphi[μ,ν] = ∂_μ∂_νφ."""
function box_coupling(bg::BGPoint, H, D1, Z, dphi, ddphi)
    gi = bg.gi; Γ = bg.Gam
    hess = [ddphi[μ, ν] - sum(Γ[λ, μ, ν] * dphi[λ] for λ in 1:4) for μ in 1:4, ν in 1:4]   # ∇_μ∇_νφ
    up = gi * hess * transpose(gi)                                                        # ∇^μ∇^νφ
    K = zeros(ComplexF64, NJ)
    for μ in 1:4, ν in 1:4
        up[μ, ν] != 0 && (K .+= up[μ, ν] .* H[μ, ν, :])
    end
    dphiup = gi * dphi
    for ν in 1:4
        K .+= dphiup[ν] .* Z[ν, :]
    end
    return K
end
