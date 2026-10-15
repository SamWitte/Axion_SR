const NO_MAIN = true
include("run_exchange.jl")
part = parse(Int, ARGS[1])
base = (a=0.0, Lextra=4, ord=8, method=:pinned, R1fac=1.1, R2fac=3.5, csfac=4.0, drho=0.02)
cases = [(0.035, (tail=50.0,)), (0.035, (tail=50.0, drho=0.015)), (0.05, (tail=50.0,)), (0.035, (tail=70.0,)),
         (0.1, (tail=50.0,)), (0.2, (tail=50.0,))]
for (al, c) in cases[part:2:end]
    println("=== a=0 alpha=$al ", c); flush(stdout)
    run_channel((3, 2, 2), (2, 1, 1), true, al; merge(base, c)...)
    flush(stdout)
end
