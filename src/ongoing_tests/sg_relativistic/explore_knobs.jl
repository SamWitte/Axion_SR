const NO_MAIN = true
include("run_exchange.jl")
part = parse(Int, ARGS[1])
base = (a=0.0, Lextra=4, ord=8, method=:pinned, R1fac=1.1, R2fac=3.5, csfac=4.0, drho=0.015)
cases = [(dkind=:d2,), (R2fac=8.0, csfac=10.0), (ord=10,), (delta=0.5,), (Lextra=7,), (nstrut=8,)]
for c in cases[part:2:end]
    println("=== a=0 alpha=0.035 ", c); flush(stdout)
    run_channel((3, 2, 2), (2, 1, 1), true, 0.035; merge(base, c)...)
    flush(stdout)
end
