const NO_MAIN = true
include("run_exchange.jl")
spin = parse(Float64, ARGS[2])
cases = [(Lextra=4, drho=0.04), (Lextra=6, drho=0.03), (Lextra=4, drho=0.025), (Lextra=8, drho=0.04)]
part = parse(Int, ARGS[1])
for c in cases[part:2:end]
    println("=== a=$spin ", c); flush(stdout)
    try
        run_channel((3, 2, 2), (2, 1, 1), true, 0.2; a=spin, c...)
    catch e
        println("FAILED: ", sprint(showerror, e))
    end
    flush(stdout)
end
