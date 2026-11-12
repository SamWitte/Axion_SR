const NO_MAIN = true
include("run_exchange.jl")
part = parse(Int, ARGS[1])
cases = [(drho=0.03,), (drho=0.025,), (drho=0.02,), (drho=0.015,)]
for c in cases[part:2:end]
    for al in (0.1, 0.05)
        println("=== a=0 alpha=$al wide window ", c); flush(stdout)
        kw = merge((a=0.0, Lextra=4, ord=8, method=:pinned, R1fac=1.1, R2fac=3.5, csfac=4.0), c)
        run_channel((3, 2, 2), (2, 1, 1), true, al; kw...)
        flush(stdout)
    end
end
