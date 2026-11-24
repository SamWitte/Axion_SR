"""
Evolution helpers for the axion-BH ODE system: boundary enforcement, spin clamping, ODE setup.
Included after Core/constants.jl.
"""

"""
    get_clamped_spin(spin_val::Float64)::Tuple{Float64, Float64}

Clamp the spin to [0, maxSpin] and return (a, rP) with rP = 1 + sqrt(1 - a^2).
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

Clamp level i in place (log u and linear u_real) to e_init <= E <= bn_list[i], the bosenova cap.
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

enforce_bosenova_boundary! for levels 1:idx_lvl.
"""
function enforce_all_boundaries!(u::Vector, u_real::Vector, idx_lvl::Int, bn_list::Vector, e_init::Float64)
    for i in 1:idx_lvl
        enforce_bosenova_boundary!(u, u_real, i, bn_list, e_init)
    end
end

"""
    check_energy_floor(u::Vector, idx_lvl::Int)::Vector{Int}

Indices of levels whose log energy u[i] is below log(SOLVER_TOLERANCES.energy_floor).
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

True if log level u[i] is within SOLVER_TOLERANCES.bosenova_threshold of log_bn.
"""
function is_near_bosenova(u::Vector, i::Int, log_bn::Float64)::Bool
    abs(u[i] - log_bn) < SOLVER_TOLERANCES.bosenova_threshold
end

"""
    initialize_solver_tolerances(non_rel::Bool, high_p::Bool)::Tuple{Float64, Float64}

Return (default_reltol, reltol_Thres), tighter for non_rel and for high_p.
"""
function initialize_solver_tolerances(non_rel::Bool, high_p::Bool)::Tuple{Float64, Float64}
    default_reltol = non_rel ? 1e-8 : 1e-7
    reltol_Thres = high_p ? 1e-7 : 1e-6
    return (default_reltol, reltol_Thres)
end

"""
    setup_state_vectors(idx_lvl::Int, aBH::Float64, M_BH::Float64, e_init::Float64, default_reltol::Float64)::Tuple{Vector, Vector}

Log-space y0 = log.([e_init (idx_lvl times); aBH; M_BH]) and reltol (default_reltol; 1e-10 for spin and mass).
"""
function setup_state_vectors(idx_lvl::Int, aBH::Real, M_BH::Real, e_init::Real, default_reltol::Real)::Tuple{Vector, Vector}
    y0 = Float64[]
    reltol = Float64[]

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
    def_spin_tol = 1e-10
    append!(reltol, def_spin_tol)
    append!(reltol, def_spin_tol)

    return (y0, reltol)
end

"""
    setup_quantum_levels_standard(Nmax::Int, fa::Float64, M_Pl::Float64, alph::Float64, aBH::Float64)::Tuple{Int, Vector, Vector, Vector}

Levels 1 <= m <= l < n <= Nmax, plus end states (2l+1, 2l, 2m) beyond Nmax, with bosenova caps e2_maxBN (n/2)^4.
"""
function setup_quantum_levels_standard(Nmax::Int, fa::Float64, M_Pl::Float64, alph::Float64, aBH::Float64)::Tuple{Int, Vector, Vector, Vector}

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

    # Truncation modes: for each (nn, l, m) with nn <= Nmax, add end state (2l+1, 2l, 2m)
    # if that state falls outside the standard levels (n_end > Nmax)
    seen_truncation_modes = Set()
    for nn in 1:Nmax, l in 1:(nn - 1), m in 1:l
        n_end = 2 * l + 1
        l_end = 2 * l
        m_end = 2 * m
        trunc_key = (n_end, l_end, m_end)
        if n_end > Nmax && !(trunc_key in seen_truncation_modes)
            idx_lvl += 1
            max_alph = aBH * m_end / (2 * (1 + sqrt(1 - aBH^2))) * 1.1
            push!(modes, (n_end, l_end, m_end, max_alph))
            push!(m_list, m_end)
            push!(bn_list, e2_maxBN * (n_end / 2)^4)
            push!(seen_truncation_modes, trunc_key)
        end
    end
    # println(idx_lvl)
    # println(m_list)
    # println(modes)
    return (idx_lvl, m_list, bn_list, modes)
end

"""
    setup_quantum_levels_spinone()::Tuple{Int, Vector, Vector, Vector}

Spinone mode: one level, idx_lvl = 1, m_list = [1], empty bn_list and modes.
"""
function setup_quantum_levels_spinone()::Tuple{Int, Vector, Vector, Vector}
    idx_lvl = 1
    m_list = [1]
    bn_list = []  # Not used in spinone mode
    modes = []     # Not used in spinone mode
    return (idx_lvl, m_list, bn_list, modes)
end
