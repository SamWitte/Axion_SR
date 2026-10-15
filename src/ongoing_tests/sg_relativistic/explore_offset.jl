const NO_MAIN = true
include("run_exchange.jl")
part = parse(Int, ARGS[1])
cases = [(frame_scale=0.0,), (R1fac=0.5, R2fac=0.8), (Lextra=6,), (R1fac=1.1, R2fac=1.45, csfac=3.0)]
for c in cases[part:2:end]
    println("=== a=0 alpha=0.05 ", c); flush(stdout)
    kw = merge((a=0.0, Lextra=4, drho=0.03, ord=8, method=:pinned), c)
    r = run_channel((3, 2, 2), (2, 1, 1), true, 0.05; kw...)
    flush(stdout)
end
