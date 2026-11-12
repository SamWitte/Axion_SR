const NO_MAIN = true
include("run_exchange.jl")
part = parse(Int, ARGS[1])
cases = [(0.025, 0.02), (0.025, 0.015), (0.035, 0.02), (0.035, 0.015)]
for (al, dr) in cases[part:2:end]
    println("=== a=0 alpha=$al drho=$dr wide window"); flush(stdout)
    for fs in (1.0, 0.0)
        run_channel((3, 2, 2), (2, 1, 1), true, al; a=0.0, Lextra=4, ord=8, method=:pinned, R1fac=1.1, R2fac=3.5, csfac=4.0, drho=dr, frame_scale=fs)
        flush(stdout)
    end
end
