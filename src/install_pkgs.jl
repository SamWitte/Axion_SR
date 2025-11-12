using Pkg

function install_packages(pkgs::Vector{String})
    # Get currently installed packages
    deps = keys(Pkg.dependencies())
    installed_pkgs = Set([string(name) for name in deps])

    for pkg in pkgs
        if pkg ∈ installed_pkgs
            println("$pkg already installed.")
        else
            println("Installing $pkg...")
            try
                Pkg.add(pkg)
            catch e
                @warn "Failed to install $pkg" exception=e
            end
        end
    end
end

install_packages(["DelimitedFiles", "PyCall", "Interpolations", "HypergeometricFunctions", "SpecialFunctions", "Dates", "Plots", "Printf", "Interpolations",
"Distributions", "Statistics", "OrdinaryDiffEq", "Random", "MCMCDiagnosticTools", "Suppressor", "KernelDensity", "StatsPlots", "Turing",
"WignerSymbols", "LinearAlgebra", "NPZ", "SpinWeightedSpheroidalHarmonics", "HypergeometricFunctions", "NLsolve", "ForwardDiff", "Glob", "ArgParse",
"QuadGK", "StatsBase", "DifferentialEquations"])

