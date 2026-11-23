# Level energies (and relativistic radial functions for n <= REL_NMAX) at a list of black-hole spins: the precompute step of
# Compute_sg_points.jl with the spin as a parameter.  Fills the gf_radial cache (same keys as production, which contain the
# spin) and writes the energies as text for the cheap spin estimator (spin_cheap.jl).
#   julia spin_levels.jl <levels file> <alphas a,b,..> <spins a,b,..> <cache dir> <out dir> <worker i> <workers N> [erg|wf|both]
# Text rows (out dir/lev_<i>of<N>.txt, appended):  n l m  spin alpha  Re(w) Im(w)  ok  seconds
#   w in units of 1/GM (compare with alpha); ok = 0 if the eigenvalue search failed or returned 0 (gf_radial then falls back to
#   the hydrogenic + fine-structure energy ergL, which is what is written in that case).
const CODE = get(ENV, "SPIN_TAB_CODE", normpath(joinpath(@__DIR__, "..", "..")))   # repository src/ (the campaign ran from a frozen copy of it)
include(joinpath(CODE, "solve_sr_rates.jl"))
include(joinpath(CODE, "Core/constants.jl"))
include(joinpath(CODE, "state_utils.jl"))
using Printf

const LEVELS, CACHE, OUT = ARGS[1], ARGS[4], ARGS[5]
const ALPHAS = parse.(Float64, split(ARGS[2], ","))
const SPINS = parse.(Float64, split(ARGS[3], ","))
const IW, NW = parse(Int, ARGS[6]), parse(Int, ARGS[7])
const WHAT = length(ARGS) >= 8 ? ARGS[8] : "both"
GF_CACHE_DIR[] = CACHE
const M_UNIT = 1.0
const NTOT_SAFE = 30000
const NPTS_BND = 2000
const REL_NMAX = 8

mkpath(OUT)
levels = [parse_state_from_args(String(strip(ln))) for ln in eachline(LEVELS) if !isempty(strip(ln)) && !startswith(strip(ln), "#")]
jobs = [(lv, al, a) for a in SPINS for lv in levels for al in ALPHAS]
key(n, l, m, a, al) = @sprintf("%d %d %d %.6g %.6g", n, l, m, a, al)
done = Set{String}()
for f in readdir(OUT; join=true)
    startswith(basename(f), "lev_") || continue
    for ln in eachline(f)
        p = split(ln)
        length(p) == 9 && push!(done, join(p[1:5], " "))
    end
end
fOUT = joinpath(OUT, "lev_$(IW)of$(NW).txt")

if WHAT in ("erg", "both")
    for (k, ((n, l, m), alpp, a)) in enumerate(jobs)
        (k - IW) % NW == 0 || continue
        key(n, l, m, a, alpp) in done && continue
        mu = alpp / (M_UNIT * GNew)
        t0 = time(); ok = 1; wr = 0.0; wi = 0.0
        try
            eG, eI = gf_cached("erg", (n, l, m, mu, M_UNIT, a, NTOT_SAFE)) do
                find_im_part(mu, M_UNIT, a, n, l, m; Ntot_force=NTOT_SAFE, for_s_rates=true, return_both=true)
            end
            wr = Float64(real(eG)); wi = Float64(real(eI))
            wr == 0 && (ok = 0)
        catch e
            ok = 0
            println("erg FAILED ", (n, l, m), " a ", a, " alpha ", alpp, ": ", first(sprint(showerror, e), 200))
        end
        ok == 0 && (wr = Float64(real(ergL(n, l, m, mu, M_UNIT, a; full=true) * GNew * M_UNIT)))
        open(fOUT, "a") do f
            println(f, key(n, l, m, a, alpp), " ", @sprintf("%.17g %.9g %d %.1f", wr, wi, ok, time() - t0))
        end
        empty!(GF_CACHE_MEM)
    end
    println("ERG_DONE ", IW, "/", NW); flush(stdout)
end

if WHAT in ("wf", "both")
    for (k, ((n, l, m), alpp, a)) in enumerate(jobs)
        (k - IW) % NW == 0 || continue
        n <= REL_NMAX || continue
        mu = alpp / (M_UNIT * GNew)
        t0 = time()
        try
            gf_cached("wf", (n, l, m, mu, M_UNIT, a, NPTS_BND, NTOT_SAFE, true)) do
                solve_radial(mu, M_UNIT, a, n, l, m; rpts=NPTS_BND, return_erg=true, Ntot_safe=NTOT_SAFE, use_heunc=true, debug=false, pre_compute_erg=nothing)
            end
        catch e
            println("wf FAILED ", (n, l, m), " a ", a, " alpha ", alpp, ": ", first(sprint(showerror, e), 200))
        end
        dt = time() - t0
        dt > 1 && (println("wf ", (n, l, m), " a ", a, " alpha ", alpp, "  ", round(dt, digits=1), " s"); flush(stdout))
        empty!(GF_CACHE_MEM)
    end
    println("WF_DONE ", IW, "/", NW); flush(stdout)
end
println("LEVELS_DONE ", IW, "/", NW)
