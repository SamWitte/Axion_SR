const NO_MAIN = true
include("run_exchange.jl")
part = parse(Int, ARGS[1])
cases = [(frame_scale=0.0, drho=0.03), (frame_scale=0.0, drho=0.025), (frame_scale=1.0, drho=0.03, R1fac=0.8, R2fac=1.2), (frame_scale=1.0, drho=0.025, R1fac=0.8, R2fac=1.2)]
for c in cases[part:2:end]
    println("=== a=0 ", c); flush(stdout)
    run_channel((3, 2, 2), (2, 1, 1), true, 0.2; a=0.0, Lextra=4, method=:pinned, c...)
    flush(stdout)
end
