using Pkg

function install_packages(pkgs::Vector{String})
    for pkg in pkgs
        if !Pkg.installed(pkg) !== nothing  # this check is deprecated; better to try/catch
            println("Installing $pkg...")
            try
                Pkg.add(pkg)
            catch e
                @warn "Failed to install $pkg" exception=e
            end
        else
            println("$pkg already installed.")
        end
    end
end

install_packages(["DelimitedFiles", "PyCall", "Interpolations", "HypergeometricFunctions", "SpecialFunctions", "Dates", "Plots", "Printf", "Interpolations",
"Distributions", "Statistics", "OrdinaryDiffEq", "Random", "MCMCDiagnosticTools", "Suppressor", "KernelDensity", "StatsPlots", "Turing", "MathLink",
"WignerSymbols", "LinearAlgebra", "NPZ", "SpinWeightedSpheroidalHarmonics", "HypergeometricFunctions", "NLsolve", "ForwardDiff", "Glob", "ArgParse",
"QuadGK", "StatsBase"])