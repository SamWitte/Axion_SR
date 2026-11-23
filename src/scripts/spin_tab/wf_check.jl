# Direct check of the cached relativistic radial functions (n <= 8) against spin: for every level, alpha and spin node the
# overlap with the hydrogenic function and with the relativistic function of the same level at a = 0.95, on a common grid.
#   julia wf_check.jl <levels file> <alphas a,b,..> <spins a,b,..> <cache dir> <out file> <worker i> <workers N>
# Rows: n l m alpha spin ok norm O_hyd O_095 Re(w)      ok = 1: entry accepted by the test of gf_radial / spin_fast.jl
#   norm = int |R|^2 r^2 dr;  O_x = |int conj(R) R_x r^2 dr| / sqrt(norms)   (1 = same shape).  Nothing is computed: a
#   missing cache entry is reported as ok = -1.
const CODE = get(ENV, "SPIN_TAB_CODE", normpath(joinpath(@__DIR__, "..", "..")))   # repository src/ (the campaign ran from a frozen copy of it)
include(joinpath(CODE, "solve_sr_rates.jl"))
include(joinpath(CODE, "Core/constants.jl"))
include(joinpath(CODE, "state_utils.jl"))
using Printf
const NPTS_BND = 2000; const NTOT_SAFE = 30000; const RPTS = 1500
trapz(y, x) = sum((y[2:end] .+ y[1:end-1]) .* diff(x)) / 2
function wf(n, l, m, mu, a, alph, lr)
    rl, r1, e = gf_cached("wf", (n, l, m, mu, 1.0, a, NPTS_BND, NTOT_SAFE, true)) do
        error("missing")
    end
    ok = !(isnan(real(e)) || isnan(imag(e)) || (abs(real(e) - alph) / alph > 1.0) || (abs(imag(e)) > 0.1))
    itp = LinearInterpolation(log10.(Float64.(rl)), ComplexF64.(r1), extrapolation_bc=Line())
    return ComplexF64.(itp(lr)), ok, Float64(real(e))
end
function main()
    levs = [parse_state_from_args(String(strip(l))) for l in eachline(ARGS[1]) if !isempty(strip(l))]
    levs = [x for x in levs if x[1] <= 8]
    alphas = parse.(Float64, split(ARGS[2], ",")); spins = parse.(Float64, split(ARGS[3], ","))
    GF_CACHE_DIR[] = ARGS[4]; iw, nw = parse(Int, ARGS[6]), parse(Int, ARGS[7])
    open(ARGS[5] * "_$(iw)of$(nw).txt", "w") do f
        for (k, (n, l, m)) in enumerate(levs)
            (k - iw) % nw == 0 || continue
            for al in alphas
                mu = al / GNew
                q = (2n) * (2n + 1); rmax = 2.5 * q / al^2
                r = 10 .^ range(log10(2.002), log10(rmax), RPTS); lr = log10.(r)
                hy = ComplexF64.(radial_bound_NR(n, l, m, mu, 1.0, r)); nh = trapz(abs2.(hy) .* r .^ 2, r)
                ref = nothing
                try; ref = wf(n, l, m, mu, 0.95, al, lr)[1]; catch; end
                for a in spins
                    try
                        R, ok, e = wf(n, l, m, mu, a, al, lr)
                        nr = trapz(abs2.(R) .* r .^ 2, r)
                        oh = abs(trapz(conj.(R) .* hy .* r .^ 2, r)) / sqrt(nr * nh)
                        o9 = ref === nothing ? NaN : abs(trapz(conj.(R) .* ref .* r .^ 2, r)) / sqrt(nr * trapz(abs2.(ref) .* r .^ 2, r))
                        @printf(f, "%d %d %d %g %g %d %.6e %.6f %.6f %.10f\n", n, l, m, al, a, ok, nr, oh, o9, e)
                    catch err
                        @printf(f, "%d %d %d %g %g -1 0 0 0 0\n", n, l, m, al, a)
                    end
                end
                length(GF_CACHE_MEM) > 200 && empty!(GF_CACHE_MEM)
            end
            flush(f)
        end
    end
    println("WF_DONE ", iw, "/", nw)
end
main()
