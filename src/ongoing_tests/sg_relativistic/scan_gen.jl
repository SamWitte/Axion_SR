# Worker for the exact self-gravity exchange of general channels (exchange_gen.jl).
# usage: julia scan_gen.jl <jobfile> <outdir> [newton_only]
# jobfile lines: S1 S2 S3 alpha [drho] [Lextra] [rmax_fac]     (states as 322 or 10-4-4; spin 0.95)
# Jobs are claimed with mkdir <outdir>/locks/<id> (many workers can share a job file); one result file per job in
# <outdir>/res/<id>.txt, one line per multipole l4 of mode 4:
#   S1 S2 S3 alpha spin drho Lextra rmax_fac dest m4 l4  Re,Im aSI  Re,Im aN  Re,Im aNR  Re,Im aGR  zmax mode_res mismatch shift seconds Lextra_dip (basis margin of dipole exchanges) version (rc2 = source cut from the wavefunctions)
include(joinpath(@__DIR__, "exchange_gen.jl"))
parse_state(s) = occursin('-', s) ? Tuple(parse.(Int, split(s, '-'))) : (parse(Int, s[1:1]), parse(Int, s[2:2]), parse(Int, s[3:3]))
function main()
    jobs = [split(l) for l in readlines(ARGS[1]) if !isempty(strip(l)) && !startswith(l, "#")]
    out = ARGS[2]; newton_only = length(ARGS) >= 3 && ARGS[3] == "newton_only"
    mkpath(joinpath(out, "locks")); mkpath(joinpath(out, "res")); mkpath(joinpath(out, "fail"))
    tstop = time() + parse(Float64, get(ENV, "SGX_WALL_SEC", "1e9"))
    spin_s = get(ENV, "SGX_SPIN", "0.95"); spin = parse(Float64, spin_s)       # test option: spin of the hole (production: 0.95)
    for j in jobs
        time() > tstop && break
        id = join(j, "_")
        isfile(joinpath(out, "res", id * ".txt")) && continue
        try mkdir(joinpath(out, "locks", id)) catch; continue end
        try write(joinpath(out, "locks", id, "owner"), get(ENV, "SLURM_JOB_ID", "none") * " " * gethostname() * "\n") catch end
        S1, S2, S3, al = j[1:4]
        dr = length(j) >= 5 ? parse(Float64, j[5]) : 0.02
        Lx = length(j) >= 6 ? parse(Int, j[6]) : 6
        rf = length(j) >= 7 ? parse(Float64, j[7]) : 1.0
        println("=== ", join(j, " ")); flush(stdout)
        try
            r = run_channel_gen(parse_state(S1), parse_state(S2), parse_state(S3), parse(Float64, al); a=spin, drho=dr, Lextra=Lx,
                                rmax_fac=rf, newton_only=newton_only, Lfull=get(ENV, "SGX_LFULL", "0") == "1")
            tmp = joinpath(out, "res", id * ".tmp." * string(getpid()))
            open(tmp, "w") do io
                for o in r.out
                    c(z) = string(real(z), " ", imag(z))
                    println(io, join([S1, S2, S3, al, spin_s, dr, Lx, rf, r.to_inf ? "Inf" : "BH", r.m4, o.l4, c(o.aSI), c(o.aN), c(o.aNR),
                                      c(o.aGR), r.zmax, r.mode_res, r.mismatch, r.shift, round(r.time, digits=1), r.Lextra_dip, r.version], " "))
                end
            end
            mv(tmp, joinpath(out, "res", id * ".txt"); force=true)
            println("  done in ", round(r.time, digits=1), " s   maxrss ", round(Sys.maxrss() / 2^30, digits=1), " GB")
        catch e
            println("FAILED: ", first(sprint(showerror, e), 400))
            open(io -> println(io, first(sprint(showerror, e), 2000)), joinpath(out, "fail", id * ".txt"), "w")
        end
        flush(stdout)
        GC.gc()                                   # return the point's arrays before the next one (nodes run short of memory)
    end
end
main()
