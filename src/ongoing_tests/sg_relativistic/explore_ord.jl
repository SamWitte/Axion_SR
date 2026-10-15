const NO_MAIN = true
include("run_exchange.jl")
part = parse(Int, ARGS[1])
cases = [(ord=8, nstrut=3, drho=0.03), (ord=8, nstrut=3, drho=0.025), (ord=6, nstrut=10, drho=0.03), (ord=6, nstrut=10, drho=0.025),
         (ord=8, nstrut=10, drho=0.03), (ord=8, nstrut=10, drho=0.02)]
for c in cases[part:2:end]
    println("=== a=0 ", c); flush(stdout)
    run_channel((3, 2, 2), (2, 1, 1), true, 0.2; a=0.0, Lextra=4, method=:pinned, R1fac=0.8, R2fac=1.2, c...)
    flush(stdout)
end
