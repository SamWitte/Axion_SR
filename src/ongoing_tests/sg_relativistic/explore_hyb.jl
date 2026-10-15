const NO_MAIN = true
include("run_exchange.jl")
spin = parse(Float64, ARGS[2])
cases = [(Lextra=4, drho=0.04), (Lextra=4, drho=0.03), (Lextra=6, drho=0.03), (Lextra=4, drho=0.025)]
part = parse(Int, ARGS[1])
for c in cases[part:2:end]
    println("=== a=$spin ", c); flush(stdout)
    try
        t0 = time()
        run_channel((3, 2, 2), (2, 1, 1), true, 0.2; a=spin, method=:pinned, c...)
        println("  time ", round(time() - t0, digits=1))
    catch e
        println("FAILED: ", sprint(showerror, e)); Base.show_backtrace(stdout, catch_backtrace())
    end
    flush(stdout)
end
