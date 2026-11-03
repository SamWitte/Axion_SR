# Accretion checks; writes data to test_plots/ (plot with plot_accretion_checks.py).
#   julia accretion_checks.jl [Nmax=3]
#   (a) Bardeen spin-up: alpha = 0.005 (SR negligible), a0 = 0.01, f_Edd = 1
#   (b) SR + accretion: M0 = 10, a0 = 0.9, alpha = 0.1, Nmax = 3, for several
#       (f_edd, track_alpha), Nmax from the command line
using DelimitedFiles, Suppressor
@suppress include(joinpath(@__DIR__, "..", "super_rad.jl"))
out = joinpath(@__DIR__, "test_plots")
Nmax = isempty(ARGS) ? 3 : parse(Int, ARGS[1])
sfx = Nmax == 3 ? "" : "_Nmax$(Nmax)"
opts = (return_all_info=true, non_rel=false, Nmax=Nmax, impose_low_cut=1e-100, eq_threshold=1e-100, abstol=1e-30)
M0 = 10.0
t, st, _, sp, mb = solve_system(0.005 / (GNew * M0), 1e18, 0.01, M0, 5e7; n_times=500, f_edd=1.0, opts...)
writedlm(joinpath(out, "bardeen$(sfx).dat"), hcat(t, sp, mb))
for (fe, ta, name) in ((0.0, false, "fedd0"), (0.0, true, "fedd0_track"), (1e-6, nothing, "fedd1e-6"),
                       (0.05, nothing, "fedd0.05"), (0.5, nothing, "fedd0.5"))
    t0 = time()
    local t, st, modes, sp, mb = solve_system(0.1 / (GNew * M0), 1e18, 0.9, M0, 1e8; n_times=2000,
                                          f_edd=fe, track_alpha=ta, opts...)
    writedlm(joinpath(out, "sr_$(name)$(sfx).dat"), hcat(t, sp, mb, hcat(st...)))
    writedlm(joinpath(out, "modes$(sfx).dat"), [m[k] for m in modes, k in 1:3])
    println(name, ": a_end=", sp[end], " M_end=", mb[end], "  (", round(time() - t0, digits=1), " s)")
end
