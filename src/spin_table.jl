# Spin dependence of the scattering rates beyond the r+ factor (NL_SPIN_TABLE, 2026-11-20).
#
# The scattering tables (rate_sve/*_LvrHc_*.dat) are computed at spin 0.95. For channels whose fourth quantum goes into
# the hole the loader and solve_system rescale them with r+(a)/r+(0.95), which is only the alpha -> 0 limit. The
# residual factor
#     R(alpha, a) = [Gamma_tab(alpha, a) / Gamma_tab(alpha, 0.95)] / rho(a),   rho = r+(a)/r+(0.95) (BH label), 1 (Inf label)
# (Gamma_tab: the tabulated quantity, i.e. without kH for BH channels with m4 != 0) is tabulated per channel for the
# self-interaction and the self-gravity amplitude separately, in rate_sve/spin_table_v1.txt or the file named by the
# environment variable NL_SPIN_TABLE (NL_SPIN_TABLE=0: no table, the r+ scaling only, as before 2026-11-23):
#     # spins 0 0.3 0.5 0.7 0.8 0.9 0.95 0.99 0.998
#     S1 S2 S3 S4 alpha flag  lnRsi_1 .. lnRsi_N  lnRsg_1 .. lnRsg_N
# (labels as in the rate tables; flag = provenance, not used here; natural logarithms). With the table
#     Gamma(a) = Gamma_SI c^2 R_SI + 2 kappa c sqrt(Gamma_SI Gamma_SG R_SI R_SG) + Gamma_SG R_SG
# (kappa from the 0.95 table), still multiplied by the r+ and kH factors of solve_system. Interpolation: monotone cubic
# (PCHIP) in a on ln R, linear in ln(alpha) between the rows of a channel, nearest row outside them; spin clamped to the
# range of the nodes; a channel without rows has R = 1.
#
# In solve_system the factor is evaluated on every RHS call as a correction to the rate used so far,
#     rate = rate_0.95 + A (R_SI - 1) + X (sqrt(R_SI R_SG) - 1) + G (R_SG - 1),
# with A, X, G the three terms above at the current mass (load_rate_coeffs(...; parts = Dict())). No ratio to the total
# is formed (it nearly cancels for some channels at f_a ~ 1e17 GeV), and a table with R = 1 leaves the run unchanged to
# the bit.

if !@isdefined(SPIN_TABLE_CACHE)
    const SPIN_TABLE_CACHE = Dict{String, Any}()
end

# The table in use. Default since 2026-11-23: rate_sve/spin_table_v1.txt (unzip rate_sve/spin_table_v1_20261123.zip there).
# NL_SPIN_TABLE=<file> selects another table; NL_SPIN_TABLE=0 (or off / none) switches the factors off (r+ scaling only,
# as before). If the default file is not there the run goes on without the factors and says so.
@isdefined(table_read) || include(joinpath(@__DIR__, "Core", "table_io.jl"))
const SPIN_TABLE_DEFAULT = joinpath(@__DIR__, "rate_sve", "spin_table_v1.txt")
if !@isdefined(SPIN_TABLE_WARNED)
    const SPIN_TABLE_WARNED = Ref(false)
end
function spin_table_path()
    p = String(strip(get(ENV, "NL_SPIN_TABLE", "")))
    lowercase(p) in ("0", "off", "none", "false") && return ""
    isempty(p) || return p
    for k in 1:3                                   # a loaded file system can report an existing file as missing
        table_exists(SPIN_TABLE_DEFAULT) && return SPIN_TABLE_DEFAULT
        k < 3 && sleep(1.0)
    end
    if !SPIN_TABLE_WARNED[]
        println("SPIN TABLE NOT FOUND: ", SPIN_TABLE_DEFAULT, " is missing (unzip rate_sve/spin_table_v1_20261123.zip into ",
                "rate_sve); the scattering rates into the hole scale with r+ only. NL_SPIN_TABLE=0 silences this.")
        SPIN_TABLE_WARNED[] = true
    end
    return ""
end

# canonical channel name, independent of the label style ("211" or "2-1-1")
function spin_table_key(s1, s2, s3, s4)
    f(s) = join(parse_state_string(String(s)), "-")
    return f(s1) * "|" * f(s2) * "|" * f(s3) * "|" * String(s4)
end
# the same from a rate key "S1_S2^S3^BH" (nothing for other key shapes)
function spin_table_key(k::AbstractString)
    p = split(k, "_", limit=2)
    length(p) == 2 || return nothing
    c = split(p[2], "^")
    (length(c) == 3 && c[3] in ("BH", "Inf")) || return nothing
    return spin_table_key(p[1], c[1], c[2], c[3])
end

# The table, read once per session: (spins, alpha of every row, V = 2 ns x rows matrix of [ln R_SI; ln R_SG],
# Dict(channel => range of its rows, sorted in alpha)).
function load_spin_table(path::AbstractString)
    get!(SPIN_TABLE_CACHE, String(path)) do
        table_exists(path) || error("NL_SPIN_TABLE: no file $(path)")
        spins = Float64[]; ns = 0
        ks = String[]; al = Float64[]; vals = Float64[]
        for ln in table_read(readlines, path)
            p = split(ln)
            isempty(p) && continue
            if startswith(p[1], "#")
                if length(p) >= 3 && p[2] == "spins"
                    spins = parse.(Float64, p[3:end]); ns = length(spins)
                end
                continue
            end
            ns == 0 && error("NL_SPIN_TABLE $(path): the '# spins ...' header must come before the rows")
            length(p) == 6 + 2ns || error("NL_SPIN_TABLE $(path): expected $(6 + 2ns) columns, got $(length(p)): $(first(ln, 80))")
            push!(ks, spin_table_key(p[1], p[2], p[3], p[4])); push!(al, parse(Float64, p[5]))
            for q in 7:(6 + 2ns)
                v = parse(Float64, p[q])
                isfinite(v) || error("NL_SPIN_TABLE $(path): non-finite value in row $(join(p[1:5], " "))")
                push!(vals, v)
            end
        end
        (ns >= 2 && issorted(spins) && allunique(spins)) || error("NL_SPIN_TABLE $(path): need a '# spins' header with >= 2 increasing nodes")
        o = sortperm(collect(zip(ks, al)))
        V = reshape(vals, 2ns, length(ks))[:, o]; ks = ks[o]; al = al[o]
        index = Dict{String, UnitRange{Int}}()
        i = 1
        while i <= length(ks)
            j = i
            while j < length(ks) && ks[j+1] == ks[i]
                j += 1
                al[j] > al[j-1] || error("NL_SPIN_TABLE $(path): repeated alpha $(al[j]) in channel $(ks[i])")
            end
            index[ks[i]] = i:j
            i = j + 1
        end
        (spins, al, V, index)
    end
end

# node slopes of the monotone cubic interpolant (Fritsch-Carlson; same formulas as scipy's PchipInterpolator)
function pchip_edge_slope(h0, h1, m0, m1)
    d = ((2h0 + h1) * m0 - h0 * m1) / (h0 + h1)
    if sign(d) != sign(m0)
        d = 0.0
    elseif sign(m0) != sign(m1) && abs(d) > 3 * abs(m0)
        d = 3 * m0
    end
    return d
end
function pchip_slopes(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    h = diff(x); m = diff(y) ./ h
    n == 2 && return [m[1], m[1]]
    d = zeros(n)
    for k in 2:n-1
        if m[k-1] * m[k] > 0
            w1 = 2h[k] + h[k-1]; w2 = h[k] + 2h[k-1]
            d[k] = (w1 + w2) / (w1 / m[k-1] + w2 / m[k])
        end
    end
    d[1] = pchip_edge_slope(h[1], h[2], m[1], m[2])
    d[n] = pchip_edge_slope(h[n-1], h[n-2], m[n-1], m[n-2])
    return d
end

# Everything the RHS kernels of solve_system need (see spin_rate). One instance per solve_system call.
mutable struct SpinFactor
    on::Bool
    spins::Vector{Float64}
    off::Vector{Int}            # per rate_cache entry: column of its first row in C (0: no row, R = 1)
    nrw::Vector{Int}            # number of rows kept (those that bracket the alpha range of the run)
    lna::Vector{Float64}        # ln(alpha) of each column of C
    C::Matrix{Float64}          # 4 ns x rows: values and slopes of ln R_SI, then of ln R_SG, at the spin nodes
    lnA::Matrix{Float64}        # ln of the SI term of each channel at the mass nodes (two equal columns without track_alpha)
    lnG::Matrix{Float64}        # ln of the SG term
    kap::Matrix{Float64}        # cross term / (2 sqrt(SI SG))
    lnMlo::Float64; dlnM::Float64; nM::Int; lnGmu::Float64; lnal0::Float64; track::Bool
    # state of the current RHS call (spin_prepare!)
    ok::Bool
    i::Int; b1::Float64; b2::Float64; b3::Float64; b4::Float64
    lnal::Float64
    j::Int; w::Float64
end
SpinFactor() = SpinFactor(false, Float64[], Int[], Int[], Float64[], zeros(0, 0), zeros(0, 0), zeros(0, 0), zeros(0, 0),
                          0.0, 1.0, 2, 0.0, 0.0, false, false, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 1, 0.0)

# parts[j][key] = (SI, cross, SG) terms of the rate at mass node j (as returned by load_rate_coeffs(...; parts)), already
# in the units of the rates the kernels use except for norm_ln[kk, j], which is added to ln SI and ln SG.
# alpha_lo, alpha_hi: range of alpha the run can reach (rows outside it are dropped, keeping one on either side).
function build_spin_factor(path, rate_keys, parts, norm_ln, alpha_lo, alpha_hi; lnMlo=0.0, dlnM=1.0, lnGmu=0.0, track=false, verbose=true)
    spins, tal, V, index = load_spin_table(path)
    ns = length(spins); nch = length(rate_keys); nM = length(parts)
    off = zeros(Int, nch); nrw = zeros(Int, nch); lna = Float64[]; cols = Vector{Float64}[]
    n_with = 0
    for (kk, k) in enumerate(rate_keys)
        ck = spin_table_key(k)
        r = ck === nothing ? nothing : get(index, ck, nothing)
        r === nothing && continue
        n_with += 1
        i1 = first(r) - 1 + something(findlast(q -> tal[q] <= alpha_lo, r), 1)
        i2 = first(r) - 1 + something(findfirst(q -> tal[q] >= alpha_hi, r), length(r))
        all(==(0.0), @view V[:, i1:i2]) && continue       # R = 1 over the whole range: nothing to do
        off[kk] = length(cols) + 1; nrw[kk] = i2 - i1 + 1
        for q in i1:i2
            push!(lna, log(tal[q]))
            push!(cols, vcat(V[1:ns, q], pchip_slopes(spins, V[1:ns, q]), V[ns+1:2ns, q], pchip_slopes(spins, V[ns+1:2ns, q])))
        end
    end
    C = isempty(cols) ? zeros(4ns, 0) : reduce(hcat, cols)
    lnA = fill(log(1e-300), nch, nM); lnG = fill(log(1e-300), nch, nM); kap = zeros(nch, nM)
    for (kk, k) in enumerate(rate_keys), j in 1:nM
        off[kk] == 0 && continue
        A, X, G = get(parts[j], k, (0.0, 0.0, 0.0))
        (isfinite(A) && isfinite(X) && isfinite(G) && abs(A + X + G) < 1e20) || continue
        A > 0 && (lnA[kk, j] = log(A) + norm_ln[kk, j])
        G > 0 && (lnG[kk, j] = log(G) + norm_ln[kk, j])
        (A > 0 && G > 0) && (kap[kk, j] = X / (2 * sqrt(A * G)))
    end
    if verbose
        println("NL_SPIN_TABLE = ", path, ": spin nodes ", spins, "; ", n_with, " of ", nch, " scattering channels have a row (",
                count(>(0), off), " with R != 1 in the alpha range of the run), ", nch - n_with, " have none (R = 1)")
    end
    return SpinFactor(true, spins, off, nrw, lna, C, lnA, lnG, kap, lnMlo, dlnM, nM, lnGmu, log(alpha_lo), track,
                      false, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 1, 0.0)
end

# Per RHS call: spin interval and Hermite weights (shared by all channels), ln(alpha) and mass bracket.
function spin_prepare!(sp::SpinFactor, a::Float64, M::Float64)
    sp.on || return
    s = sp.spins; n = length(s)
    sp.ok = isfinite(a)
    sp.ok || return
    ac = clamp(a, s[1], s[n])
    i = clamp(searchsortedlast(s, ac), 1, n - 1)
    h = s[i+1] - s[i]; t = (ac - s[i]) / h
    sp.i = i
    sp.b1 = (1 + 2t) * (1 - t)^2; sp.b2 = h * t * (1 - t)^2; sp.b3 = t^2 * (3 - 2t); sp.b4 = h * t^2 * (t - 1)
    if sp.track
        x = (log(M) - sp.lnMlo) / sp.dlnM
        x = isnan(x) ? 0.0 : clamp(x, 0.0, Float64(sp.nM - 1))       # as mass_bracket in solve_system
        j = clamp(floor(Int, x) + 1, 1, sp.nM - 1)
        sp.j = j; sp.w = clamp(x - (j - 1), 0.0, 1.0)
        sp.lnal = sp.lnGmu + sp.lnMlo + sp.dlnM * x
    else
        sp.j = 1; sp.w = 0.0; sp.lnal = sp.lnal0
    end
    return
end

# (ln R_SI, ln R_SG) of channel kk at the spin and alpha of the last spin_prepare!
@inline function spin_lnR(sp::SpinFactor, kk::Int)
    o = sp.off[kk]
    (o == 0 || !sp.ok) && return 0.0, 0.0
    nr = sp.nrw[kk]; la = sp.lnal
    r = o; wa = 0.0
    if nr > 1 && la > sp.lna[o]
        if la >= sp.lna[o+nr-1]
            r = o + nr - 1
        else
            while sp.lna[r+1] < la
                r += 1
            end
            wa = (la - sp.lna[r]) / (sp.lna[r+1] - sp.lna[r])
        end
    end
    C = sp.C; ns = length(sp.spins); i = sp.i
    @inbounds begin
        lsi = sp.b1 * C[i, r] + sp.b2 * C[ns+i, r] + sp.b3 * C[i+1, r] + sp.b4 * C[ns+i+1, r]
        lsg = sp.b1 * C[2ns+i, r] + sp.b2 * C[3ns+i, r] + sp.b3 * C[2ns+i+1, r] + sp.b4 * C[3ns+i+1, r]
        if wa != 0.0
            q = r + 1
            lsi2 = sp.b1 * C[i, q] + sp.b2 * C[ns+i, q] + sp.b3 * C[i+1, q] + sp.b4 * C[ns+i+1, q]
            lsg2 = sp.b1 * C[2ns+i, q] + sp.b2 * C[3ns+i, q] + sp.b3 * C[2ns+i+1, q] + sp.b4 * C[3ns+i+1, q]
            lsi = (1 - wa) * lsi + wa * lsi2
            lsg = (1 - wa) * lsg + wa * lsg2
        end
    end
    return lsi, lsg
end

# Rate of channel kk with the spin factors applied; base = the rate the kernels used so far (0.95 tables).
@inline function spin_rate(sp::SpinFactor, kk::Int, base::Float64)
    (base == 0.0 || sp.off[kk] == 0) && return base
    lsi, lsg = spin_lnR(sp, kk)
    (lsi == 0.0 && lsg == 0.0) && return base
    j = sp.j; w = sp.w
    @inbounds begin
        lA = (1 - w) * sp.lnA[kk, j] + w * sp.lnA[kk, j+1]
        lG = (1 - w) * sp.lnG[kk, j] + w * sp.lnG[kk, j+1]
        ka = (1 - w) * sp.kap[kk, j] + w * sp.kap[kk, j+1]
    end
    r = base + exp(lA) * (exp(lsi) - 1) + 2 * ka * exp(0.5 * (lA + lG)) * (exp(0.5 * (lsi + lsg)) - 1) + exp(lG) * (exp(lsg) - 1)
    return (r * base < 0) ? 0.0 : r      # the corrected rate keeps the sign of the tabulated one
end
