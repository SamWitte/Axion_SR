"""
    Evolution helper functions for axion-black hole ODE system.

Provides utilities for boundary condition enforcement, spin clamping, and ODE setup.
These functions are included into files that have already included Constants.jl
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

# Details
Identifies levels where quantum coherence is effectively lost (underflow regime)
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
