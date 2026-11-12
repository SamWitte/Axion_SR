const NO_MAIN = true
include("run_exchange.jl")
cases = [(zeta=1.0, wfac=30.0, csfac=1.5), (zeta=0.5, wfac=30.0, csfac=1.5), (zeta=2.0, wfac=30.0, csfac=1.5),
         (zeta=1.0, wfac=60.0, csfac=1.5), (zeta=1.0, wfac=15.0, csfac=1.5), (zeta=1.0, wfac=30.0, csfac=3.0)]
part = parse(Int, ARGS[1])
for c in cases[part:2:end]
    println("=== ", c); flush(stdout)
    run_channel((3, 2, 2), (2, 1, 1), true, 0.2; Lextra=4, drho=0.04, c...)
    flush(stdout)
end
