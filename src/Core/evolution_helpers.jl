"""
    Evolution helper functions for axion-black hole ODE system.

Provides utilities for boundary condition enforcement, spin clamping, and ODE setup.
These functions are included into files that have already included Core/constants.jl
"""

"""
    get_clamped_spin(spin_val::Float64)::Tuple{Float64, Float64}

Clamp black hole spin to valid range and compute Kerr metric factor rP.

# Arguments
- `spin_val::Float64`: Current black hole spin (0 ≤ a ≤ 1)

# Returns
- `(clamped_spin::Float64, rP::Float64)`: Clamped spin and rP = 1 + sqrt(1 - a²)

# Details
Enforces spin bounds: 0 ≤ a ≤ 0.998 (maxSpin)
Special handling at boundaries for Kerr geometry stability.
"""
function get_clamped_spin(spin_val::Float64)::Tuple{Float64, Float64}
    if spin_val > maxSpin
        return (maxSpin, 1.0 + sqrt(1.0 - maxSpin^2))
    elseif spin_val < 0.0
        return (0.0, 2.0)
    else
        return (spin_val, 1.0 + sqrt(1.0 - spin_val^2))
    end
end

"""
    enforce_bosenova_boundary!(u::Vector, u_real::Vector, i::Int, bn_list::Vector, e_init::Float64)

Enforce bosenova (binding energy) boundary conditions on quantum level i.

# Arguments
- `u::Vector`: Log-space state vector (modified in place)
- `u_real::Vector`: Real-space state vector (modified in place)
- `i::Int`: Index of quantum level to constrain
- `bn_list::Vector`: List of bosenova binding energies
- `e_init::Float64`: Minimum energy threshold

# Details
Ensures binding energy stays within physically valid bounds:
- e_init ≤ E ≤ bn_list[i] (bosenova maximum)
Uses log-space comparison with SOLVER_TOLERANCES.bosenova_threshold tolerance
"""
function enforce_bosenova_boundary!(u::Vector, u_real::Vector, i::Int, bn_list::Vector, e_init::Float64)
    # Check if below minimum energy
    if u_real[i] < e_init
        u_real[i] = e_init
        u[i] = log(e_init)
    end

    # Check if at or above bosenova threshold
    if (abs(u[i] - log(bn_list[i])) < SOLVER_TOLERANCES.bosenova_threshold) || (u[i] > log(bn_list[i]))
        u[i] = log(bn_list[i])
        u_real[i] = bn_list[i]
    end
end

"""
    enforce_all_boundaries!(u::Vector, u_real::Vector, idx_lvl::Int, bn_list::Vector, e_init::Float64)

Apply boundary condition enforcement to all quantum levels.

# Arguments
- `u::Vector`: Log-space state vector (modified in place)
- `u_real::Vector`: Real-space state vector (modified in place)
- `idx_lvl::Int`: Number of active quantum levels
- `bn_list::Vector`: List of bosenova binding energies
- `e_init::Float64`: Minimum energy threshold

# Details
Calls enforce_bosenova_boundary! for each level i ∈ [1, idx_lvl]
"""
function enforce_all_boundaries!(u::Vector, u_real::Vector, idx_lvl::Int, bn_list::Vector, e_init::Float64)
    for i in 1:idx_lvl
        enforce_bosenova_boundary!(u, u_real, i, bn_list, e_init)
    end
end

"""
    check_energy_floor(u::Vector, idx_lvl::Int)::Vector{Int}

Identify quantum levels that have fallen below energy floor.

# Arguments
- `u::Vector`: Log-space state vector
- `idx_lvl::Int`: Number of active quantum levels

# Returns
- `Vector{Int}`: Indices of levels below energy floor (e < 1e-75)

"""
function check_energy_floor(u::Vector, idx_lvl::Int)::Vector{Int}
    below_floor = Int[]
    for i in 1:idx_lvl
        if u[i] < log(SOLVER_TOLERANCES.energy_floor)
            push!(below_floor, i)
        end
    end
    return below_floor
end

"""
    is_near_bosenova(u::Vector, i::Int, log_bn::Float64)::Bool

Check if quantum level is near bosenova threshold (within tolerance).

# Arguments
- `u::Vector`: Log-space state vector
- `i::Int`: Index of level to check
- `log_bn::Float64`: Log of bosenova binding energy for this level

# Returns
- `Bool`: true if level is within SOLVER_TOLERANCES.bosenova_threshold of bosenova
"""
function is_near_bosenova(u::Vector, i::Int, log_bn::Float64)::Bool
    abs(u[i] - log_bn) < SOLVER_TOLERANCES.bosenova_threshold
end

"""
    initialize_solver_tolerances(non_rel::Bool, high_p::Bool)::Tuple{Float64, Float64}

Compute default and threshold tolerances for ODE solver based on physical regime.

# Arguments
- `non_rel::Bool`: Whether using non-relativistic approximation
- `high_p::Bool`: Whether using high-precision mode

# Returns
- `(default_reltol::Float64, reltol_Thres::Float64)`: Default and threshold relative tolerances

# Details
Non-relativistic mode uses tighter tolerances (1e-5 vs 1e-3).
High-precision mode uses much tighter threshold (1e-3 vs 0.1).
"""
function initialize_solver_tolerances(non_rel::Bool, high_p::Bool)::Tuple{Float64, Float64}
    default_reltol = non_rel ? 1e-5 : 1e-3
    reltol_Thres = high_p ? 1e-3 : 0.1
    return (default_reltol, reltol_Thres)
end

"""
    setup_state_vectors(idx_lvl::Int, aBH::Float64, M_BH::Float64, e_init::Float64, default_reltol::Float64)::Tuple{Vector, Vector}

Initialize ODE state vector y0 and relative tolerance array reltol.

# Arguments
- `idx_lvl::Int`: Number of active quantum levels
- `aBH::Float64`: Black hole spin
- `M_BH::Float64`: Black hole mass
- `e_init::Float64`: Initial energy threshold
- `default_reltol::Float64`: Default relative tolerance for energy levels

# Returns
- `(y0::Vector, reltol::Vector)`: Log-space initial state and per-component tolerances

# Details
Creates state vector [log(e1), ..., log(e_idx_lvl), log(aBH), log(M_BH)].
Applies default_reltol to energy components and def_spin_tol=1e-3 to spin/mass.
"""
function setup_state_vectors(idx_lvl::Int, aBH::Real, M_BH::Real, e_init::Real, default_reltol::Real)::Tuple{Vector, Vector}
    y0 = []
    reltol = []

    # Energy levels initialization
    for i in 1:idx_lvl
        append!(y0, e_init)
        append!(reltol, default_reltol)
    end

    # Spin and mass initialization
    append!(y0, aBH)
    append!(y0, M_BH)
    y0 = log.(y0)

    # Tighter tolerance for spin and mass
    def_spin_tol = 1e-3
    append!(reltol, def_spin_tol)
    append!(reltol, def_spin_tol)

    return (y0, reltol)
end

"""
    setup_quantum_levels_standard(Nmax::Int, fa::Float64, M_Pl::Float64, alph::Float64, aBH::Float64)::Tuple{Int, Vector, Vector, Vector}

Setup quantum levels for standard (non-spinone) mode with multiple levels.

# Arguments
- `Nmax::Int`: Maximum principal quantum number (3-8)
- `fa::Float64`: Axion decay constant (1/GeV)
- `M_Pl::Float64`: Planck mass
- `alph::Float64`: Dimensionless axion-gravity coupling
- `aBH::Float64`: Black hole spin

# Returns
- `(idx_lvl::Int, m_list::Vector, bn_list::Vector, modes::Vector)`: Number of levels, m values, bosenova energies, mode tuples

# Details
Creates all (n,l,m) levels with 1 ≤ m ≤ l, 1 ≤ l ≤ n-1, 1 ≤ n ≤ Nmax.
Also adds truncation modes for high-m states.
Bosenova binding energy: e_bn = e2_maxBN * (n/2)^4 with e2_maxBN = 1024π(fa/M_Pl)²/(9α³)
"""
function setup_quantum_levels_standard(Nmax::Int, fa::Float64, M_Pl::Float64, alph::Float64, aBH::Float64)::Tuple{Int, Vector, Vector, Vector}
    if !(3 <= Nmax <= 8)
        error("Nmax must be between 3 and 8, got $Nmax")
    end

    idx_lvl = 0
    m_list = []
    bn_list = []
    modes = []
    truncation_modes = []

    # Compute maximum bosenova energy scale
    e2_maxBN = 1024 * pi * (fa / M_Pl)^2 / (9 * alph^3)

    # Standard (n,l,m) levels
    for nn in 1:Nmax, l in 1:(nn - 1), m in 1:l
        idx_lvl += 1
        max_alph = aBH * m / (2 * (1 + sqrt(1 - aBH^2))) * 1.1  # Growth estimate with 10% safety margin
        push!(modes, (nn, l, m, max_alph))
        push!(m_list, m)
        push!(bn_list, e2_maxBN * (nn / 2)^4)
    end

    # Truncation modes for high-m states
    for nn in 1:Nmax, l in 1:(nn - 1)
        m = l
        m_new = 2 * m
        if (m_new >= Nmax) && !any(mode -> mode[1:3] == (m_new + 1, m_new, m_new), truncation_modes)
            idx_lvl += 1
            max_alph = aBH * m_new / (2 * (1 + sqrt(1 - aBH^2))) * 1.1
            push!(modes, (m_new + 1, m_new, m_new, max_alph))
            push!(truncation_modes, (m_new + 1, m_new, m_new, max_alph))
            push!(m_list, m_new)
            push!(bn_list, e2_maxBN * ((m_new + 1) / 2)^4)
        end
    end

    return (idx_lvl, m_list, bn_list, modes)
end

"""
    setup_quantum_levels_spinone()::Tuple{Int, Vector, Vector, Vector}

Setup quantum levels for spinone mode (single level only).

# Returns
- `(idx_lvl::Int, m_list::Vector, bn_list::Vector, modes::Vector)`: idx_lvl=1, m=[1], empty lists for others

# Details
Spinone mode tracks only a single quantum level with m=1.
No bosenova energy tracking needed (handled differently).
"""
function setup_quantum_levels_spinone()::Tuple{Int, Vector, Vector, Vector}
    idx_lvl = 1
    m_list = [1]
    bn_list = []  # Not used in spinone mode
    modes = []     # Not used in spinone mode
    return (idx_lvl, m_list, bn_list, modes)
end
