# usage: julia scan.jl <jobfile> <part> <nparts> <outfile>
# jobfile lines: S1 S3 dest alpha spin drho Lextra
const NO_MAIN = true
include("run_exchange.jl")
jobs = [split(l) for l in readlines(ARGS[1]) if !isempty(strip(l)) && !startswith(l, "#")]
part = parse(Int, ARGS[2]); np = parse(Int, ARGS[3]); out = ARGS[4]
st(s) = (parse(Int, s[1]), parse(Int, s[2]), parse(Int, s[3]))
for (k, j) in enumerate(jobs)
    (k - 1) % np + 1 == part || continue
    S1, S3, dest, al, sp, dr, Lx = j
    println("=== ", join(j, " ")); flush(stdout)
    t0 = time()
    try
        r = run_channel(st(S1), st(S3), dest == "Inf", parse(Float64, al); a=parse(Float64, sp), drho=parse(Float64, dr),
                        Lextra=parse(Int, Lx))
        open(out, "a") do io
            println(io, join(j, " "), " ", r.l4, " ", real(r.gr_over_n), " ", imag(r.gr_over_n), " ", real(r.nr_over_n), " ",
                    imag(r.nr_over_n), " ", r.C_N, " ", real(r.si_over_n), " ", imag(r.si_over_n), " ", round(time() - t0, digits=1))
        end
    catch e
        println("FAILED: ", sprint(showerror, e)); Base.show_backtrace(stdout, catch_backtrace()); println()
    end
    flush(stdout)
end
