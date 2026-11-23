# Klein-Gordon charge Q of every level of the campaign at every alpha of the grid, with the radial function and the
# frequency that gf_radial uses for it (Kerr solution from the cache for n <= 8, hydrogenic function above; flat
# normalisation int |R|^2 r^2 dr = 1), computed with kg_norm_ratio of solve_sr_rates.jl.  Used by allW/merge_allW.py to
# convert the stored rate points (computed with the flat normalisation) to the Klein-Gordon normalisation:
#     rate -> rate / (Q_1 Q_2 Q_3).
#   julia kg_levels.jl <levels file> <alpha list> <cache dir> <out dir> <worker i> <N workers>
const SRC = "/data/dust/user/samwitte/Axion_SR_sg/src"
include(joinpath(SRC, "solve_sr_rates.jl")); include(joinpath(SRC, "Core/constants.jl")); include(joinpath(SRC, "state_utils.jl"))
using Printf
const LEV, ALS, CACHE, OUT = ARGS[1], parse.(Float64, split(ARGS[2], ",")), ARGS[3], ARGS[4]
const IW, NW = parse(Int, ARGS[5]), parse(Int, ARGS[6])
GF_CACHE_DIR[] = CACHE
const A_SPIN = 0.95; const M_UNIT = 1.0; const NTOT_SAFE = 30000; const NPTS_BND = 2000; const REL_NMAX = 8
levels = [String(strip(l)) for l in eachline(LEV) if !isempty(strip(l))]
jobs = [(s, al) for s in levels for al in ALS]; mkpath(OUT)
open(joinpath(OUT, "part_$(IW).txt"), "a") do io
    for (k, (s, alpp)) in enumerate(jobs)
        (k - IW) % NW == 0 || continue
        n, l, m = parse_state_from_args(s); a = A_SPIN; M = M_UNIT; mu = alpp / (M * GNew); alph = alpp
        rp = 1 + sqrt(1 - a^2)
        rmax = n < 10 ? 2.0^(2.0 * n - 2 * (1 + n)) * gamma(2 + 2 * n) / alph^2 / factorial(2 * n - 1) * 10.0 : 2.0^(2.0 * n - 2 * (1 + n)) * 4 * n^2 / alph^2 * 10.0
        rlist = 10 .^ range(log10(rp * (1.0 + 1e-3)), log10(rmax), 3000)
        kind = "nr"; q = NaN; wre = NaN
        try
            eG = gf_cached("erg", (n, l, m, mu, M, a, NTOT_SAFE)) do
                find_im_part(mu, M, a, n, l, m; Ntot_force=NTOT_SAFE, for_s_rates=true, return_both=true)
            end
            ergG = eG[1] == 0 ? ergL(n, l, m, mu, M, a; full=true) * GNew * M : eG[1]
            rf = nothing; erg = ergG
            if n <= REL_NMAX
                rl, r1, e1 = gf_cached("wf", (n, l, m, mu, M, a, NPTS_BND, NTOT_SAFE, true)) do
                    solve_radial(mu, M, a, n, l, m; rpts=NPTS_BND, return_erg=true, Ntot_safe=NTOT_SAFE, use_heunc=true, debug=false, pre_compute_erg=nothing)
                end
                if !(isnan(real(e1)) || isnan(imag(e1)) || (abs(real(e1) - alph) / alph > 1.0) || (abs(imag(e1)) > 0.1))
                    itp = LinearInterpolation(log10.(rl), r1, extrapolation_bc=Line()); rf = itp(log10.(rlist)); erg = e1; kind = "kerr"
                end
            end
            rf === nothing && (rf = radial_bound_NR(n, l, m, mu, M, rlist))
            wre = Float64(real(erg)); q = kg_norm_ratio(rf, rlist, l, m, a, alph, wre)
        catch e
            println("FAILED ", s, " alpha ", alpp, ": ", first(sprint(showerror, e), 160))
        end
        @printf(io, "%s %.10g %.10f %s %.10f\n", s, alpp, q, kind, wre / alph); flush(io)
    end
end
println("KG_LEVELS_DONE worker $IW")
