# Exact one-graviton exchange for general channels  N1 x N2 -> N3 x {BH, Inf}  (N1 and N2 may differ, any n).
#
# The source for mode 4 is the symmetrised sum of the two exchanges, exactly as in the Newtonian pipeline
# (gf_radial(...; self_grav=true): 1/2 [psi1 Phi_23 + psi2 Phi_13]):
#     J_GR = 1/2 [ S4(h[T_23] acting on phi1) + S4(h[T_13] acting on phi2) ] / (-4 alpha^4)
# with T_b3 the transition stress tensor of b -> 3 (frequency w_b - w3, azimuthal number m_b - m3), h[T] the
# Lorenz-gauge solution on Kerr (solver.jl; hole-frame "pinned" solver for |m_b - m3| = 1) and S4 = -delta box[h] phi_a.
# For N1 = N2 this is the calculation of run_exchange.jl (one exchange).
#
# Everything is projected on mode 4 with the Kerr Green's function of the production rates (gf_amplitudes), for every
# multipole l4 the production code sums (BH: l4 = lmin, lmin + 2, ... <= BHlmax; Inf: the single leading l4).
# Output per l4: the amplitudes of the contact (SI), Newtonian (fixed hole), Newtonian + recoil and exact sources.
#
# Differences from run_exchange.jl:
#  * bound states from the radial ODE (OdeMode below) instead of the Leaver series.  The series sum_k b_k u^k cancels
#    catastrophically for states with radial nodes (terms / sum ~ (r/GM)^(n-l-1): 1e16 for 10-4-4), so in Float64 it
#    is only usable for the nodeless states l = n - 1.  The eigenfrequency still comes from the continued fraction.
#  * assembly and the mode-4 source are threaded over the radial index (JULIA_NUM_THREADS).
#  * the angular basis of each exchange is cut at min(l_b + l3, l_a + max l4) + Lextra: only multipoles L <= l_a + l4
#    of the potential of b -> 3 can take mode a to mode 4.
#
# usage: see scan_gen.jl
const NO_MAIN = true
include("run_exchange.jl")

# ---------------------------------------------------------------- bound states from the radial ODE
# Ingoing-regular radial function f(r) = R_BL(r) e^{i w r* - i m r#} of the mode (real frequency w):
#     Delta f'' + (Delta' - 2 i K) f' - (i K' + lamT + mu^2 r^2) f = 0,   K = (r^2 + a^2) w - a m,
#     lamT = Lam + a^2 w^2 - 2 a m w,   (L^2 - a^2 (w^2 - mu^2) cos^2) S = Lam S.
# Inner part: f integrated outwards from the horizon (regular solution, grows through the barrier: stable).
# Outer part: the real BL function R integrated inwards from kappa r = tailmax (decaying solution: stable).
# The two are joined at r_m = n^2/alpha^2, inside the classically allowed region; the mismatch of the logarithmic
# derivatives there measures the error of the eigenfrequency.
struct OdeMode
    w::Float64; m::Int; a::Float64; alph::Float64; lamT::Float64
    r0::Float64; rm::Float64; rout::Float64
    ein::Vector{Float64}; sin::Vector{Any}          # segment edges (increasing r) and solutions of the inner part
    eout::Vector{Float64}; sout::Vector{Any}        # segment edges (decreasing r) and solutions of the outer part
    s::ComplexF64
    tay::NTuple{4,ComplexF64}
    mismatch::Float64
end
const ODE_MODES = IdDict{Any,OdeMode}()

# Integrate u' = f(u, r) through the segments edges[k] -> edges[k+1].  The absolute tolerance of each segment is tol
# times the size of the solution at its start: the solution changes by many orders of magnitude and goes through
# zero at the radial nodes, where a purely relative tolerance stalls the step-size control.
function solve_segments(f!, u0, edges, kref; tol=1e-12)
    sols = Any[]; u = copy(u0)
    for k in 1:length(edges)-1
        kl = kref(edges[k])
        sc = abs(u[1]) + abs(u[2]) / kl
        sol = solve(ODEProblem(f!, u, (edges[k], edges[k+1])), Vern9(); reltol=tol, abstol=[tol * sc, tol * sc * kl], maxiters=10^7)
        push!(sols, sol); u = copy(sol.u[end])
    end
    return sols
end

function build_ode_mode(n, m, alph, a, w, Lam; tailmax=80.0, tol=1e-12)
    rp = rplus(a); lamT = real(Lam) + a^2 * w^2 - 2a * m * w
    Dl(r) = r^2 - 2r + a^2
    Kf(r) = (r^2 + a^2) * w - a * m
    Cf(r) = 2im * r * w + lamT + alph^2 * r^2
    Bf(r) = 2 * (r - 1) - 2im * Kf(r)
    # Taylor coefficients at the horizon (f = 1): successive derivatives of the ODE at Delta = 0
    Dp = 2rp - 2; Bp = 2 - 4im * rp * w; Cp = 2im * w + 2 * alph^2 * rp; Bpp = -4im * w; Cpp = 2 * alph^2
    t1 = Cf(rp) / Bf(rp)
    t2 = (Cp + (Cf(rp) - Bp) * t1) / (Dp + Bf(rp))
    t3 = (Cpp + (2Cp - Bpp) * t1 + (Cf(rp) - 2Bp - 2) * t2) / (2Dp + Bf(rp))
    eps = 2e-5; r0 = rp + eps
    u0 = ComplexF64[1 + t1 * eps + t2 * eps^2 / 2 + t3 * eps^3 / 6, t1 + t2 * eps + t3 * eps^2 / 2]
    rm = max(n^2 / alph^2, 3.0)
    fin!(du, u, p, r) = (du[1] = u[2]; du[2] = (Cf(r) * u[1] - Bf(r) * u[2]) / Dl(r); nothing)
    nin = ceil(Int, log((rm - rp) / eps) / log(1.2))
    ein = rp .+ eps .* ((rm - rp) / eps) .^ (collect(0:nin) ./ nin)
    sin = solve_segments(fin!, u0, ein, r -> max(alph, 1 / (r - rp + 1)); tol=tol)
    kap = sqrt(alph^2 - w^2)
    rout = max(tailmax / kap, 3rm)
    fout!(du, u, p, r) = (D = Dl(r); du[1] = u[2]; du[2] = -(2r - 2) / D * u[2] - (Kf(r)^2 / D^2 - (lamT + alph^2 * r^2) / D) * u[1]; nothing)
    nout = ceil(Int, (rout - rm) * kap * 2)
    eout = collect(range(rout, rm, length=nout + 1))
    kloc(r) = (D = Dl(r); sqrt(abs(Kf(r)^2 / D^2 - (lamT + alph^2 * r^2) / D)) + 1 / r)
    sout = solve_segments(fout!, [1.0, -kap], eout, kloc; tol=tol)
    # join at rm: compare (R, R'/k) of the two solutions in BL form
    D = Dl(rm); kl = kloc(rm)
    ph = exp(im * (w * rstar(a, rm) - m * rsharp(a, rm)))
    fi, fi1 = sin[end].u[end]
    vi = [fi / ph, (fi1 - im * Kf(rm) / D * fi) / ph / kl]
    Ro, Ro1 = sout[end].u[end]
    vo = ComplexF64[Ro, Ro1 / kl]
    s = sum(conj.(vi) .* vo) / sum(abs2.(vi))
    mismatch = sqrt(sum(abs2.(s .* vi .- vo)) / sum(abs2.(vo)))
    return OdeMode(w, m, a, alph, lamT, r0, rm, rout, ein, sin, eout, sout, s, (1.0 + 0im, t1, t2, t3), mismatch)
end

function ode_radial_f(om::OdeMode, r::Real)
    a = om.a; w = om.w; m = om.m; rp = rplus(a)
    r > om.rout && return (0.0im, 0.0im, 0.0im)
    if r < om.r0
        x = r - rp; t = om.tay
        return (om.s * (t[1] + t[2] * x + t[3] * x^2 / 2 + t[4] * x^3 / 6), om.s * (t[2] + t[3] * x + t[4] * x^2 / 2), om.s * (t[3] + t[4] * x))
    end
    D = r^2 - 2r + a^2; K = (r^2 + a^2) * w - a * m
    if r <= om.rm
        k = clamp(searchsortedlast(om.ein, r), 1, length(om.sin))
        u = om.sin[k](r); f = om.s * u[1]; f1 = om.s * u[2]
    else
        k = clamp(searchsortedlast(om.eout, r; rev=true), 1, length(om.sout))
        u = om.sout[k](r); ph = exp(im * (w * rstar(a, r) - m * rsharp(a, r)))
        f = u[1] * ph; f1 = (u[2] + im * K / D * u[1]) * ph
    end
    f2 = ((2im * r * w + om.lamT + om.alph^2 * r^2) * f - (2 * (r - 1) - 2im * K) * f1) / D
    return (f, f1, f2)
end

# series evaluation of modes.jl (kept for modes without an ODE table and for complex r)
function radial_f_series(md::KerrMode, r)
    a = md.a; rp, rm = rplus(a), rminus(a); L = md.L; w = md.w; m = md.m
    p = im * L.sig + L.chi - 1 - im * (2rm * w - m * a) / (rp - rm)
    u = (r - rp) / (r - rm); du = (rp - rm) / (r - rm)^2; d2u = -2 * (rp - rm) / (r - rm)^3
    S0 = 0.0im; S1 = 0.0im; S2 = 0.0im
    for k in length(md.b):-1:1
        S2 = S2 * u + 2 * S1
        S1 = S1 * u + S0
        S0 = S0 * u + md.b[k]
    end
    lp = p / (r - rm) + (L.q + im * w)
    P = exp(p * log(r - rm) + (L.q + im * w) * r)
    P1 = P * lp; P2 = P * (lp^2 - p / (r - rm)^2)
    return md.norm * P * S0, md.norm * (P1 * S0 + P * S1 * du), md.norm * (P2 * S0 + 2 * P1 * S1 * du + P * (S2 * du^2 + S1 * d2u))
end

# radial_f of modes.jl, redirected to the ODE table when the mode has one
function radial_f(md::KerrMode, r)
    om = get(ODE_MODES, md, nothing)
    (om === nothing || !(r isa Real)) && return radial_f_series(md, r)
    f, f1, f2 = ode_radial_f(om, r)
    return md.norm * f, md.norm * f1, md.norm * f2
end

"""Kerr mode with the eigenfrequency from the continued fraction (as kerr_mode) and radial functions from the ODE."""
function kerr_mode_ode(n, l, m, alph, a, w0; tol=1e-13, maxit=60)
    Kt = cf_terms(alph, w0, a)
    function resid(w)
        Lam, _ = spheroidal_coeffs(l, m, a^2 * (w^2 - alph^2))
        return leaver_cf(leaver_coeffs(w, alph, a, m, Lam), Kt)[2]
    end
    w1 = complex(w0); w2 = w1 * (1 + 1e-7) + 1e-12im
    f1 = resid(w1); f2 = resid(w2)
    for it in 1:maxit
        w3 = w2 - f2 * (w2 - w1) / (f2 - f1)
        w1, f1 = w2, f2
        w2 = w3; f2 = resid(w2)
        abs(w2 - w1) < tol * abs(w2) && break
    end
    w = w2
    Lam, sph = spheroidal_coeffs(l, m, a^2 * (w^2 - alph^2))
    om = build_ode_mode(n, m, alph, a, real(w), Lam)
    md0 = KerrMode(n, l, m, alph, a, w, Lam, sph, leaver_coeffs(w, alph, a, m, Lam), ComplexF64[1.0], 1.0)
    ODE_MODES[md0] = om
    rr = 10 .^ range(log10(rplus(a) * (1 + 1e-6)), log10(om.rout), 6000)
    fv = [abs2(radial_f(md0, r)[1]) for r in rr]
    nm = sum(0.5 .* (fv[1:end-1] .* rr[1:end-1] .^ 2 .+ fv[2:end] .* rr[2:end] .^ 2) .* diff(rr))
    md = KerrMode(n, l, m, alph, a, w, Lam, sph, md0.L, md0.b, 1 / sqrt(nm))
    ODE_MODES[md] = om
    return md, abs(f2), om.mismatch
end

const MODE_CACHE = Dict{Any,Any}()
function get_mode(st, alph, a)
    get!(MODE_CACHE, (st, alph, a)) do
        n, l, m = st
        g = find_im_part(alph / GNew, 1.0, a, n, l, m; Ntot_force=5000, for_s_rates=true, return_both=true)
        w0 = ComplexF64(Float64(g[1]) + im * Float64(g[2]))
        md, res, mis = kerr_mode_ode(n, l, m, alph, a, w0)
        # the refinement must stay on the level found by find_im_part (level spacing ~ alpha^3 / n^3)
        shift = abs(real(md.w) - real(w0)) / (alph^3 / n^3)
        (md, res, mis, shift)
    end
end

kbound(md) = real(sqrt(md.alph^2 - real(md.w)^2 + 0im))

# ---------------------------------------------------------------- threaded assembly and sources
# Rows of radial point i written straight into the CSC arrays of the transpose A^T (column (i-1) n + u of A^T is row
# (i-1) n + u of A: dense n x n blocks for the radial points of the finite-difference stencil).  Each call writes its
# own range, so the calls can run in threads; no intermediate triplets (memory = the matrix itself).
function fill_rows_T!(rowval, nzval, colptr, cols, i, a, w, m, ab::AngBasis, rg::RadGrid, Hj)
    n = NC(ab); N = length(rg.r)
    if i == N
        for u in 1:n
            p = colptr[(i-1)*n+u]; rowval[p] = (i - 1) * n + u; nzval[p] = 1.0
        end
        return nothing
    end
    A0, A1, A2 = angular_blocks(a, rg.r[i], w, m, ab, Hj; kap=r -> 0.0)
    rp1 = rg.dr[i]; rp2 = rg.d2r[i]
    Bρ = A1 ./ rp1 .- A2 .* (rp2 / rp1^3)
    Bρρ = A2 ./ rp1^2
    d1 = rg.D1[i, :]; d2 = rg.D2[i, :]
    blk = zeros(ComplexF64, n, n)
    for (kc, jn) in enumerate(cols)
        blk .= d1[jn] .* Bρ .+ d2[jn] .* Bρρ
        jn == i && (blk .+= A0)
        for u in 1:n
            p0 = colptr[(i-1)*n+u] + (kc - 1) * n - 1
            for v in 1:n
                rowval[p0+v] = (jn - 1) * n + v; nzval[p0+v] = blk[u, v]
            end
        end
    end
    return nothing
end

"""Transpose of the global operator of assemble() (solver.jl), threaded over the radial index."""
function assemble_T(a, w, m, ab::AngBasis, rg::RadGrid)
    n = NC(ab); N = length(rg.r)
    Hj = jet_arrays(w, m)
    colsets = [i == N ? [N] : sort(union(findnz(rg.D1[i, :])[1], findnz(rg.D2[i, :])[1], [i])) for i in 1:N]
    colptr = Vector{Int}(undef, N * n + 1); colptr[1] = 1
    for i in 1:N, u in 1:n
        c = (i - 1) * n + u
        colptr[c+1] = colptr[c] + (i == N ? 1 : length(colsets[i]) * n)
    end
    nz = colptr[end] - 1
    rowval = Vector{Int}(undef, nz); nzval = Vector{ComplexF64}(undef, nz)
    Threads.@threads for i in 1:N
        fill_rows_T!(rowval, nzval, colptr, colsets[i], i, a, w, m, ab, rg, Hj)
    end
    return SparseMatrixCSC(N * n, N * n, colptr, rowval, nzval)
end

function minus_dG_col(es, hfun, R2, Hj, i)
    ab = es.ab; r = es.rg.r[i]
    (imag(r) != 0 || real(r) >= R2) && return zeros(ComplexF64, NC(ab))
    r = real(r)
    Tc = zeros(ComplexF64, 10, length(ab.T.th))
    for (j, th) in enumerate(ab.T.th)
        bg = bg_point(es.a, r, th)
        D1, D2 = cov_derivs(bg, Hj...)
        _, dG, _ = einstein_jets(bg, D1, D2)
        J = analytic_jets(hfun, r, th)
        for (cp, (μ, ν)) in enumerate(PAIRS)
            Tc[cp, j] = -sum(dG[μ, ν, :] .* J)
        end
    end
    return project_tensor(ab, Tc)
end

function minus_dG_source_mt(es, hfun, R2)
    n = NC(es.ab); N = length(es.rg.r)
    Hj = jet_arrays(es.wg, es.mg)
    s = zeros(ComplexF64, n, N)
    Threads.@threads for i in 1:N
        s[:, i] = minus_dG_col(es, hfun, R2, Hj, i)
    end
    return s
end

blas_serial() = LinearAlgebra.BLAS.set_num_threads(1)
blas_parallel() = LinearAlgebra.BLAS.set_num_threads(max(2, Threads.nthreads()))

function solve_h_mt(es::ExchangeSetup; verbose=true)
    t0 = time()
    b = rhs_T31(es)
    blas_serial(); At = assemble_T(es.a, es.wg, es.mg, es.ab, es.rg); blas_parallel()
    t1 = time()
    nz = nnz(At)
    F = lu(At); At = nothing
    X = reshape(transpose(F) \ vec(b), NC(es.ab), length(es.rg.r))
    F = nothing; GC.gc()
    verbose && println("  assemble $(round(t1-t0,digits=1))s  solve $(round(time()-t1,digits=1))s  N=$(length(es.rg.r)) NC=$(NC(es.ab)) nnz=$nz")
    return X
end

# solve_h_pinned of pinned.jl (Kerr: odd = :frame) with the threaded assembly; the bordered system is transposed too:
# [A Z; P^T 0]^T = [A^T P; Z^T 0]
function solve_h_pinned_mt(es; R1, R2, verbose=true, nstrut=3)
    t0 = time()
    n = NC(es.ab); N = length(es.rg.r)
    b = rhs_T31(es)
    blas_serial(); At = assemble_T(es.a, es.wg, es.mg, es.ab, es.rg)
    pd, pt, zd, zt = pins_and_struts(es; nstrut=nstrut)
    sc = maximum(abs, At.nzval)
    pdn = pd .* (sc / maximum(abs, pd)); ptn = pt .* sc
    odd = es.a != 0
    P = odd ? [pdn ptn] : reshape(pdn, :, 1)
    Z = odd ? [zd zt] : reshape(zd, :, 1)
    np_ = size(P, 2)
    kinds = odd ? [:d, :t] : [:d]
    AbT = [At sparse(P); sparse(transpose(Z)) spzeros(ComplexF64, np_, np_)]
    At = nothing
    blas_parallel(); F = transpose(lu(AbT)); AbT = nothing
    t1 = time()
    fields = [windowed(k, es.mg, R1, R2) for k in kinds]
    xs = []
    x0 = F \ [vec(b); zeros(np_)]
    for f in fields
        blas_serial(); sk = minus_dG_source_mt(es, f, R2); blas_parallel()
        ck = zeros(ComplexF64, n * N); ck[1:n] = analytic_coeffs(es, f, 1)
        push!(xs, F \ [vec(sk); -(transpose(P) * ck)])
    end
    F = nothing; GC.gc()
    nk = length(kinds)
    tau(x) = x[n*N+1:n*N+nk]
    M = reduce(hcat, [tau(x) for x in xs])
    acoef = -(M \ tau(x0))
    # Test switch SGX_HOLD_FIXED=1: do not add the frame fields, i.e. leave the struts in place.  The hole is then held at
    # rest (no displacement, no tilt) by an external force instead of recoiling freely; the result must reproduce the
    # Newtonian amplitude WITHOUT recoil at small alpha.  Off by default (the physical, strut-free solution).
    get(ENV, "SGX_HOLD_FIXED", "0") == "1" && (acoef = zero(acoef))
    Xv = x0[1:n*N]; tres = x0[n*N+1:end]
    for k in 1:nk
        Xv = Xv .+ acoef[k] .* xs[k][1:n*N]
        tres = tres .+ acoef[k] .* xs[k][n*N+1:end]
    end
    hA = (i, j, r, th) -> sum(acoef[k] * fields[k](i, j, r, th) for k in 1:nk)
    verbose && println("  pinned: assemble+factor $(round(t1-t0,digits=1))s, frame+solves $(round(time()-t1,digits=1))s  N=$N NC=$n  ",
                       "residual struts / tau0 = ", round.(abs.(tres[1:nk]) ./ abs.(tau(x0)), sigdigits=2))
    return reshape(Xv, n, N), hA
end

"""Outermost radius at which r^2 g(r) exceeds tol times its maximum (g sampled on a log grid from rlo to rhi)."""
function last_above(g, rlo, rhi; tol=1e-10, npts=3000)
    rr = 10 .^ range(log10(rlo), log10(rhi), npts)
    v = [g(r) * r^2 for r in rr]
    k = findlast(x -> x > tol * maximum(v), v)
    return rr[min(k + 1, npts)]
end
mode_rout(md) = (om = get(ODE_MODES, md, nothing); om === nothing ? 80.0 / kbound(md) : om.rout)

"""Setup for the exchange b -> 3 acting on mode a.  Returns (es, Sa) with Sa = (S, dS, d2S) of mode a on the nodes."""
function setup_gen(md_b, md3, md_a; Lmax, drho=0.02, ord=8, zeta=1.0, delta=1.0, src_tol=1e-15, csfac=4.0, wfac=30.0, wid=1.5, tail=50.0, cut_tol=1e-10)
    a = md_b.a; alph = md_b.alph
    wg = real(md_b.w) - real(md3.w); mg = md_b.m - md3.m
    ab = AngBasis(mg, Lmax)
    # T_b3 is negligible beyond -log(src_tol)/(k_b + k3); the mode-4 source h grad grad phi_a must be kept until phi_a
    # itself is negligible (k_a r = tail), see setup() in exchange.jl
    rcut = max(-log(src_tol) / (kbound(md_b) + kbound(md3)), tail / kbound(md_a))
    # The two estimates above assume pure exponentials.  High-n levels oscillate out to ~2 n^2/alpha^2 before they decay,
    # so for two high-n levels b, 3 and a low-n level a they end where psi_b psi_3* is still large (322 x 15-4-4 -> 12-6-6:
    # cut at 230 a0, source extends to ~300 a0).  Keep everything up to the last radius where r^2 |R_b R_3| (source) or
    # r^2 |R_a| (acted-on level) exceeds cut_tol of its maximum.
    rcut = max(rcut, last_above(r -> abs(radial_f(md_b, r)[1]) * abs(radial_f(md3, r)[1]), rplus(a) * 1.001, min(mode_rout(md_b), mode_rout(md3)); tol=cut_tol),
               last_above(r -> abs(radial_f(md_a, r)[1]), rplus(a) * 1.001, mode_rout(md_a); tol=cut_tol))
    rcs = csfac * rcut
    rend = rcs * 3 + wfac / abs(wg)
    rg = RadGrid(a, wg; delta=delta, rcs=rcs, rend=rend, drho=drho, zeta=zeta, ord=ord, wid=wid)
    ireal = findlast(x -> imag(x) == 0 && real(x) <= rcut, rg.r)
    nth = length(ab.T.th)
    onnodes(md) = spheroidal_on_nodes(md.sph, SWSHTable(md.m, maximum(keys(md.sph)), nth; smax=0))
    es = ExchangeSetup(md_b, md3, a, alph, wg, mg, ab, rg, onnodes(md_b), onnodes(md3), ireal, delta)
    return es, onnodes(md_a)
end

function mode4_row(es::ExchangeSetup, X, fa, Sa, wa, ma, Hj, hA, i)
    ab = es.ab; rg = es.rg
    S, dS, d2S = Sa
    r = real(rg.r[i])
    row = zeros(ComplexF64, length(ab.T.th))
    zmax = 0.0; dmax = 0.0
    for (j, th) in enumerate(ab.T.th)
        bg = bg_point(es.a, r, th)
        D1, D2 = cov_derivs(bg, Hj...)
        _, _, Z = einstein_jets(bg, D1, D2)
        J = h_jets(ab, rg, X, i, j)
        if hA === nothing
            zmax = max(zmax, maximum(abs(sum(Z[ν, :] .* J)) for ν in 1:4))
            dmax = max(dmax, maximum(abs(sum(D1[ρ, μ, ν, :] .* J)) for ρ in 1:4, μ in 1:4, ν in 1:4))
        else
            J = J .+ analytic_jets(hA, r, th)
        end
        d1, dd = scalar_jets(fa..., S[j], dS[j], d2S[j], wa, ma)
        K = box_coupling(bg, Hj[1], D1, Z, d1, dd)
        row[j] = sum(K .* J)
    end
    return row, zmax / max(dmax, 1e-300)
end

"""Mode-4 source S4[i, j] = -delta box[h] phi_a on the real radial points (h = numerical X plus optional analytic hA)."""
function mode4_source_gen(es::ExchangeSetup, X, md_a, Sa; hA=nothing)
    wa = real(md_a.w); ma = md_a.m
    Hj = jet_arrays(es.wg, es.mg)
    S4 = zeros(ComplexF64, es.ireal, length(es.ab.T.th))
    Zrel = zeros(es.ireal)
    fas = [radial_f(md_a, real(es.rg.r[i])) for i in 1:es.ireal]          # serial: the ODE interpolants are not thread-safe
    Threads.@threads for i in 1:es.ireal
        row, z = mode4_row(es, X, fas[i], Sa, wa, ma, Hj, hA, i)
        S4[i, :] = row; Zrel[i] = z
    end
    return S4, Zrel
end

"""Exact exchange b -> 3 acting on a (mode 4 up to l4max): returns (es, S4, zmax)."""
function exchange_piece(md_b, md3, md_a, l4max; Lextra=6, Lextra_dip=10, Lfull=false, R1fac=1.1, R2fac=3.5, nstrut=3, verbose=true, setup_kw...)
    # dipole exchanges (hole-frame solver) need a larger basis: the recoil cancels most of the direct L = 1 term, and the
    # remainder is converged only for Lextra >= 10 (211 x 544 -> 655 BH at alpha = 0.1: A_GR/A_NR = 3.05+4.77i, 1.006+0.118i,
    # 1.011+0.005i, 1.012+0.002i for Lextra = 6, 8, 10, 12; 211 x 211 -> 322: 1.016504, 1.016462, 1.016460)
    Lx = abs(md_b.m - md3.m) == 1 ? max(Lextra, Lextra_dip) : Lextra
    Lmax = (Lfull ? md_b.l + md3.l : min(md_b.l + md3.l, md_a.l + l4max)) + Lx
    es, Sa = setup_gen(md_b, md3, md_a; Lmax=Lmax, setup_kw...)
    rc = real(es.rg.r[es.ireal])
    if abs(es.mg) == 1
        X, hA = solve_h_pinned_mt(es; R1=R1fac * rc, R2=R2fac * rc, verbose=verbose, nstrut=nstrut)
        blas_serial(); S4, _ = mode4_source_gen(es, X, md_a, Sa; hA=hA); blas_parallel()
        zmax = 0.0
    else
        X = solve_h_mt(es; verbose=verbose)
        blas_serial(); S4, Zrel = mode4_source_gen(es, X, md_a, Sa); blas_parallel()
        wts = [maximum(abs.(S4[i, :])) for i in 1:es.ireal]
        zmax = maximum(Zrel[wts .> 1e-3 * maximum(wts)])
    end
    return es, S4, zmax
end

# ---------------------------------------------------------------- Newtonian and contact sources (BL, same modes)
"""
Radial sources on the fine BL grid r, projected on conj(S4) with the Sigma weight, for mode 4 = (l4, m4, w4):
  SI : R1 R2 R3* <S1 S2 S3* S4* Sigma>
  N  : 1/2 [ R1 Phi_23 + R2 Phi_13 ]          (Newtonian, fixed hole)
  NR : same with the BH recoil in the dipole (L = 1) part
"""
function newtonian_sources_gen(md1, md2, md3, l4, m4, w4, r, Rs; nth=96, Lextra=6)
    a = md1.a; alph = md1.alph
    _, sph4 = spheroidal_coeffs(l4, m4, a^2 * (w4^2 - alph^2))
    onn(md) = spheroidal_on_nodes(md.sph, SWSHTable(md.m, maximum(keys(md.sph)), nth; smax=0))[1]
    T4 = SWSHTable(m4, maximum(keys(sph4)), nth; smax=0)
    S1 = onn(md1); S2 = onn(md2); S3 = onn(md3); S4 = spheroidal_on_nodes(sph4, T4)[1]
    x = T4.x; w = T4.w
    R1, R2, R3 = Rs
    ang(v) = 2pi * sum(w .* v)
    A0 = ang(S1 .* S2 .* conj.(S3) .* conj.(S4)); A2 = ang(S1 .* S2 .* conj.(S3) .* conj.(S4) .* x .^ 2)
    T_SI = R1 .* R2 .* conj.(R3) .* (A0 .* r .^ 2 .+ A2 * a^2)
    # Test switch SGX_NEWTON_SHIFT = d (units of GM; off = 0): solve the Poisson equation of the Newtonian term in the
    # radial coordinate x = r_BL - d instead of r_BL (kernel and volume element), everything else unchanged.  d = 1 is
    # the harmonic (also, to this order, the isotropic) radial coordinate of the hole, in which the post-Newtonian
    # expansion is written; the difference to d = 0 is of relative order GM/r ~ alpha^2, i.e. a 1PN-order statement
    # about WHICH Newtonian term the 1PN corrections of ongoing_tests/sg_vs_si/sg_1pn*.jl are corrections to.
    xs = r .- parse(Float64, get(ENV, "SGX_NEWTON_SHIFT", "0"))
    # Test switches SGX_NEWTON_WX = wx, SGX_NEWTON_WY = wy (off = 0): weight the Newtonian term (fixed hole) with the
    # potential of the hole, (1 + wx GM/x) at the point where the potential acts and (1 + wy GM/x) in the transition
    # density that sources it.  The complete 1PN source of the t/u-channel exchange without recoil, derived from the
    # Klein-Gordon equation in the 1PN harmonic metric (note, Section 3.7), is the Newtonian one in the harmonic radius
    # with wx = 4, wy = 6, times 1 + 2 (w_a w_4 + w_b w_3 - 2 mu^2)/mu^2, plus retardation and gravitomagnetism.
    wxf = 1 .+ parse(Float64, get(ENV, "SGX_NEWTON_WX", "0")) ./ xs
    wyf = 1 .+ parse(Float64, get(ENV, "SGX_NEWTON_WY", "0")) ./ xs
    # Recoil term at 1PN (same derivation, in the frame of the hole): the fictitious potential a_BH . x couples with
    # the weight (1 + 5/2 GM/x) instead of (1 + 4 GM/x) -> SGX_NEWTON_WRX = 2.5 (default: the same weight as the direct
    # term), and the acceleration of the hole gets the retardation term of the potential,
    # a_BH ~ int rho (1 + w_g^2 x^2 / 2) dx -> SGX_NEWTON_RECRET = 1.  (The gravitomagnetic induction term of a_BH is
    # computed with ongoing_tests/sg_vs_si/sg_1pn_recoil.jl.)
    wrf = haskey(ENV, "SGX_NEWTON_WRX") ? 1 .+ parse(Float64, ENV["SGX_NEWTON_WRX"]) ./ xs : wxf
    recret = get(ENV, "SGX_NEWTON_RECRET", "0") == "1"
    function piece(mdb, Sb, Rb, Sa, Ra)                      # potential of b -> 3 acting on a
        M_ = mdb.m - md3.m
        Lmax = mdb.l + md3.l + Lextra
        TL = SWSHTable(M_, Lmax, nth; smax=0)
        rho = Rb .* conj.(R3) .* wyf
        TN = zeros(ComplexF64, length(r)); TNR = zeros(ComplexF64, length(r))
        for L in abs(M_):Lmax
            cL = ang(Sb .* conj.(S3) .* TL.Y[(0, L)])
            abs(cL) < 1e-14 && continue
            ci = cumtrap(rho .* xs .^ (L + 2), r); co = revtrap(rho .* xs .^ (1 - L), r)
            Q = (xs .^ (-(L + 1)) .* ci .+ xs .^ L .* co) .* wxf
            QR = Q
            if L == 1
                co1 = co[1]
                recret && (co1 += 0.5 * (real(mdb.w) - real(md3.w))^2 * revtrap(rho .* xs .^ 2, r)[1])
                QR = Q .- xs .* co1 .* wrf
            end
            BL0 = ang(Sa .* TL.Y[(0, L)] .* conj.(S4)); BL2 = ang(Sa .* TL.Y[(0, L)] .* conj.(S4) .* x .^ 2)
            TN .+= (4pi / (2L + 1)) * cL .* Q .* (BL0 .* r .^ 2 .+ BL2 * a^2)
            TNR .+= (4pi / (2L + 1)) * cL .* QR .* (BL0 .* r .^ 2 .+ BL2 * a^2)
        end
        return Ra .* TN, Ra .* TNR
    end
    N23, NR23 = piece(md2, S2, R2, S1, R1)
    N13, NR13 = md1 === md2 ? (N23, NR23) : piece(md1, S1, R1, S2, R2)
    return T_SI, 0.5 .* (N23 .+ N13), 0.5 .* (NR23 .+ NR13), sph4
end

"""Multipoles of mode 4 summed by the production code (gf_radial with BHlmax)."""
function l4_list(st1, st2, st3, to_inf; BHlmax=2)
    l1, m1 = st1[2], st1[3]; l2, m2 = st2[2], st2[3]; l3, m3 = st3[2], st3[3]
    m4 = m1 + m2 - m3
    if to_inf
        l = l1 + l2 - l3
        while !(l >= abs(m4) && l > 0)
            l += 2
        end
        return [l], m4
    end
    l = max(abs(l1 + l2 - l3), abs(m4))
    return collect(l:2:max(BHlmax, l)), m4
end

"""
One channel at one alpha.  newton_only = true skips the graviton solves (Newtonian and contact amplitudes only).
rmax_fac scales the outer radius of the mode-4 integration (test of the truncation for Inf channels).
"""
function run_channel_gen(st1, st2, st3, alph; a=0.95, Lextra=6, drho=0.02, rpts=10000, verbose=true, newton_only=false,
                         rmax_fac=1.0, alpha_h=0.12, BHlmax=2, Lfull=false, kw...)
    t0 = time()
    same = st1 == st2
    md1, r1, q1, s1 = get_mode(st1, alph, a); md2, r2, q2, s2 = same ? (md1, r1, q1, s1) : get_mode(st2, alph, a)
    md3, r3, q3, s3 = get_mode(st3, alph, a)
    verbose && @printf("modes: w = %.12f %.12f %.12f  cf res %.1e %.1e %.1e  join mismatch %.1e %.1e %.1e  level shift %.1e %.1e %.1e\n",
                       real(md1.w), real(md2.w), real(md3.w), r1, r2, r3, q1, q2, q3, s1, s2, s3)
    w4 = real(md1.w) + real(md2.w) - real(md3.w)
    to_inf = w4 > alph
    ls, m4 = l4_list(st1, st2, st3, to_inf; BHlmax=BHlmax)

    # fine BL grid for the Green's function (as run_exchange.jl / gf_radial)
    rp = rplus(a); eps_fac = 1e-3; maxN = max(st1[1], st2[1], st3[1])
    rbnd = 0.25 * (2maxN + 1) * (2maxN) / alph^2
    rmax = (to_inf ? 30.0 / sqrt(w4^2 - alph^2) : rbnd * 10.0) * rmax_fac
    rlist = collect(10 .^ range(log10(rp * (1.0 + eps_fac)), log10(rmax), rpts))
    rmax_1 = 0.25 * 7 * 6 / 0.03^2 * 7.0
    h_mve = max(0.2 / (rbnd * 7.0 / rmax_1) * 10.0, rp * alph / alpha_h)      # step of the mode-4 integration <= alpha_h/alpha
    # test option SGX_HYDRO_NMIN = n0: hydrogenic radial functions for levels with n >= n0 in the Newtonian and contact
    # sources (what the production tables do above n = 8); off by default
    hyd0 = parse(Int, get(ENV, "SGX_HYDRO_NMIN", "0"))
    rad(st, md) = (hyd0 > 0 && st[1] >= hyd0) ? ComplexF64.(Float64.(radial_bound_NR(st[1], st[2], st[3], alph / GNew, 1.0, rlist))) :
                                                 ComplexF64[R_BL(md, ri) for ri in rlist]
    Rb1 = rad(st1, md1)
    Rs = (Rb1, same ? Rb1 : rad(st2, md2), rad(st3, md3))

    pieces = Any[]; zmax = 0.0
    if !newton_only
        es, S4, z = exchange_piece(md2, md3, md1, maximum(ls); Lextra=Lextra, Lfull=Lfull, drho=drho, verbose=verbose, kw...)      # 2 -> 3, acts on 1
        push!(pieces, (es, S4)); zmax = max(zmax, z)
        if !same
            es, S4, z = exchange_piece(md1, md3, md2, maximum(ls); Lextra=Lextra, Lfull=Lfull, drho=drho, verbose=verbose, kw...)  # 1 -> 3, acts on 2
            push!(pieces, (es, S4)); zmax = max(zmax, z)
        end
    end
    out = []
    for l4 in ls
        T_SI, T_N, T_NR, sph4 = newtonian_sources_gen(md1, md2, md3, l4, m4, w4, rlist, Rs; Lextra=Lextra)
        srcs = [T_SI, T_N, T_NR]
        if !newton_only
            T_GR = zeros(ComplexF64, length(rlist))
            for (es, S4) in pieces
                rs, Tg = gr_radial_source(es, S4, l4, m4, w4, sph4)
                T_GR .+= (same ? 1.0 : 0.5) .* interp_rho(es, rs, Tg, rlist)
            end
            push!(srcs, T_GR)
        end
        amps, _, _ = gf_amplitudes(complex(w4), a, alph, l4, m4, rlist, srcs, rmax, to_inf; h_mve=h_mve, eps_fac=eps_fac)
        aGR = newton_only ? complex(NaN, NaN) : amps[4]
        push!(out, (; l4, aSI=amps[1], aN=amps[2], aNR=amps[3], aGR))
        verbose && @printf("  alpha=%.4f l4=%d m4=%d %s: A_GR/A_NR = %.6f%+.6fi   A_NR/A_N = %.6f%+.6fi   C_NR = %.5e   |A_SI|^2 = %.4e\n",
                           alph, l4, m4, to_inf ? "Inf" : "BH", real(aGR / amps[3]), imag(aGR / amps[3]), real(amps[3] / amps[2]),
                           imag(amps[3] / amps[2]), 64 * alph^4 * abs2(amps[3]) / abs2(amps[1]), abs2(amps[1]))
    end
    pieces = nothing; GC.gc()
    return (; alph, a, to_inf, m4, w4, zmax, out, time=time() - t0, mode_res=max(r1, r2, r3), mismatch=max(q1, q2, q3), shift=max(s1, s2, s3), Lextra_dip=10, version="rc2")
end
