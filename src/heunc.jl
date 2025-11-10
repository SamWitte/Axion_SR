# heun_switch_solver.jl
using DifferentialEquations, LinearAlgebra, Logging

# ---------- Helpers: 4th-order Taylor initializer for u (ComplexF64) ----------
function taylor_init_u4_f64(q::ComplexF64, α::ComplexF64, γ::ComplexF64,
                            δ::ComplexF64, ε::ComplexF64, z0::ComplexF64, λ::ComplexF64)
    ONE = ComplexF64(1.0,0.0)
    TWO = ComplexF64(2.0,0.0)

    # y-derivatives at 0 (a0..a4), then transform to u-derivs b0..b4
    a0 = ONE
    a1 = -q/γ
    p1 = γ + δ - ε
    p2 = ε
    a2 = ((p1 - q) * a1 + α * a0) / (ONE + γ)
    denom3 = ONE + γ/2
    num3 = (ONE + p1 - q/2) * a2 + (p2 + α) * a1
    a3 = num3 / denom3
    denom4 = (3.0 + γ) / 6.0
    num4 = (ONE + p1/2 - q/6) * a3 + (p2 + α/2) * a2
    a4 = num4 / denom4

    # u-derivs: invert binomial relation with λ
    b0 = a0
    b1 = a1 - λ*b0
    b2 = a2 - (λ^2)*b0 - 2*λ*b1
    b3 = a3 - (λ^3)*b0 - 3*(λ^2)*b1 - 3*λ*b2
    b4 = a4 - (λ^4)*b0 - 4*(λ^3)*b1 - 6*(λ^2)*b2 - 4*λ*b3

    z1 = z0
    z2 = z1*z1
    z3 = z2*z1
    z4 = z3*z1

    u_z0 = b0 + b1*z1 + b2*(z2)/2 + b3*(z3)/6 + b4*(z4)/24
    uprime_z0 = b1 + b2*z1 + b3*(z2)/2 + b4*(z3)/6

    return ComplexF64[u_z0, uprime_z0]
end

# ---------- Transformed ODE: f_u_inplace and Jacobian (ComplexF64) ----------
# R(z) for transformed equation
R_of_f64(z, α, γ, δ, ε, q) = α*z - q + z*(z-1)*(ε^2/4) - (ε/2)*(γ*(z-1) + δ*z)

function f_u_inplace_f64!(out::Vector{ComplexF64}, uvec::Vector{ComplexF64}, z::ComplexF64,
                          α, γ, δ, ε, q)
    out[1] = uvec[2]
    denom = z*(z-1)
    denom_abs = abs(denom)
    tiny = 1e-24
    tiny_z = 1e-12
    if denom_abs < tiny
        if abs(z) < tiny_z
            R0 = R_of_f64(0.0+0.0im, α, γ, δ, ε, q)
            R1 = α - (ε^2/4) - (ε/2)*(γ + δ)
            u1 = -q/γ + ε/2
            u2 = ((γ + δ + R0)*u1 + R1*(1.0+0.0im)) / (1.0+γ)
            out[2] = u2
        else
            A = γ*(z-1) + δ*z
            B = R_of_f64(z, α, γ, δ, ε, q)
            reg = ComplexF64(denom_abs + tiny, 0.0)
            out[2] = (-A*uvec[2] - B*uvec[1]) / (denom + reg)
        end
    else
        A = γ*(z-1) + δ*z
        B = R_of_f64(z, α, γ, δ, ε, q)
        out[2] = (-A*uvec[2] - B*uvec[1]) / denom
    end
    return out
end

function jac_u_f64!(J::Matrix{ComplexF64}, u::Vector{ComplexF64}, z::ComplexF64, α, γ, δ, ε, q)
    denom = z*(z-1)
    denom_abs = abs(denom)
    tiny = 1e-24
    tiny_z = 1e-12
    if denom_abs < tiny
        if abs(z) < tiny_z
            jf11 = ComplexF64(0.0,0.0); jf12 = ComplexF64(0.0,0.0)
        else
            A = γ*(z-1) + δ*z
            B = R_of_f64(z, α, γ, δ, ε, q)
            reg = ComplexF64(denom_abs + tiny, 0.0)
            jf11 = -B / (denom + reg)
            jf12 = -A / (denom + reg)
        end
    else
        A = γ*(z-1) + δ*z
        B = R_of_f64(z, α, γ, δ, ε, q)
        jf11 = -B / denom
        jf12 = -A / denom
    end
    # ∂f/∂u matrix (2x2): [0 1; jf11 jf12]
    J[1,1] = ComplexF64(0.0,0.0); J[1,2] = ComplexF64(1.0,0.0)
    J[2,1] = jf11; J[2,2] = jf12
    return nothing
end

# ---------- Standard (y) ODE: f_y_inplace and Jacobian (ComplexF64) ----------
# original ODE: z(z-1) y'' + (γ(z-1)+δ z + z(z-1) ε) y' + (α z - q) y = 0
function f_y_inplace_f64!(out::Vector{ComplexF64}, yvec::Vector{ComplexF64}, z::ComplexF64,
                           γ, δ, ε, α, q)
    out[1] = yvec[2]
    denom = z*(z-1)
    denom_abs = abs(denom)
    tiny = 1e-24
    tiny_z = 1e-12
    if denom_abs < tiny
        if abs(z) < tiny_z
            # not expected: near 0; but return safe value
            out[2] = 0.0 + 0.0im
        else
            A = γ*(z-1) + δ*z + z*(z-1)*ε
            B = α*z - q
            reg = ComplexF64(denom_abs + tiny, 0.0)
            out[2] = (-A*yvec[2] - B*yvec[1]) / (denom + reg)
        end
    else
        A = γ*(z-1) + δ*z + z*(z-1)*ε
        B = α*z - q
        out[2] = (-A*yvec[2] - B*yvec[1]) / denom
    end
    return out
end

function jac_y_f64!(J::Matrix{ComplexF64}, y::Vector{ComplexF64}, z::ComplexF64, γ, δ, ε, α, q)
    denom = z*(z-1)
    denom_abs = abs(denom)
    tiny = 1e-24
    tiny_z = 1e-12
    if denom_abs < tiny
        if abs(z) < tiny_z
            jf11 = ComplexF64(0.0,0.0); jf12 = ComplexF64(0.0,0.0)
        else
            A = γ*(z-1) + δ*z + z*(z-1)*ε
            B = α*z - q
            reg = ComplexF64(denom_abs + tiny, 0.0)
            jf11 = -B / (denom + reg)
            jf12 = -A / (denom + reg)
        end
    else
        A = γ*(z-1) + δ*z + z*(z-1)*ε
        B = α*z - q
        jf11 = -B / denom
        jf12 = -A / denom
    end
    # Jacobian of first-order system du/dz = [y'; y''] is [0 1; jf11 jf12]
    J[1,1] = ComplexF64(0.0,0.0); J[1,2] = ComplexF64(1.0,0.0)
    J[2,1] = jf11; J[2,2] = jf12
    return nothing
end

# ---------- Integrators for segments (generic wrappers) ----------
# integrate a segment for a provided in-place RHS and jacobian
function integrate_segment_generic(z_start::ComplexF64, z_end::ComplexF64, u_start::Vector{ComplexF64},
                                   rhs_inplace!, jac_inplace!, params...;
                                   reltol=1e-12, abstol=1e-14)
    dz = z_end - z_start
    function rhs_dt!(du, u, p, t)
        zt = z_start + t*dz
        tmp = similar(du)
        rhs_inplace!(tmp, u, zt, params...)
        du[1] = tmp[1] * dz
        du[2] = tmp[2] * dz
        return nothing
    end
    function jac_dt!(J, u, p, t)
        zt = z_start + t*dz
        # J_f = ∂f/∂u (2x2)
        Jtmp = similar(zeros(ComplexF64,2,2))
        jac_inplace!(Jtmp, u, zt, params...)
        # multiply by dz to get Jacobian for du/dt
        J[1,1] = Jtmp[1,1] * dz; J[1,2] = Jtmp[1,2] * dz
        J[2,1] = Jtmp[2,1] * dz; J[2,2] = Jtmp[2,2] * dz
        return nothing
    end
    of = ODEFunction(rhs_dt!, jac = jac_dt!)
    prob = ODEProblem(of, u_start, (0.0, 1.0))
    sol = nothing
    try
        sol = solve(prob, RadauIIA5(autodiff=false), reltol=reltol, abstol=abstol)
    catch err
        sol = solve(prob, Rodas5(autodiff=false), reltol=reltol, abstol=abstol)
    end
    return sol(1.0)
end

# ---------- High-level routine: get u (renormalized) at target via segments ----------
# returns (u_norm, M) where actual u = M * u_norm
function integrate_u_with_renorm(q, α, γ, δ, ε, z_target::ComplexF64;
                                nseg::Int=200, imag_eps=0.0, z0_abs=1e-9, reltol=1e-12, abstol=1e-14,
                                renorm::Bool=true)

    # Set lambda
    λ = -ε/ComplexF64(2.0,0.0)

    # build complex end with optional tiny imag offset
    zend = ComplexF64(real(z_target), imag_eps)
    # directional small start near 0
    if zend == 0.0+0.0im
        zstart = ComplexF64(z0_abs, imag_eps)
    else
        amp = abs(zend)
        zstart = (amp == 0.0) ? ComplexF64(z0_abs, imag_eps) : (zend/ComplexF64(amp,0.0)) * ComplexF64(z0_abs, imag_eps)
    end

    # make nodes
    nodes = [ zstart + (k/(nseg))*(zend - zstart) for k in 0:nseg ]

    # initial u via 4th-order Taylor at zstart
    u = taylor_init_u4_f64(q, α, γ, δ, ε, zstart, λ)
    M = ComplexF64(1.0,0.0)   # cumulative multiplier

    for k in 1:nseg
        za = nodes[k]      # nodes array indexed from 1
        zb = nodes[k+1]
        # integrate segment
        u_end = integrate_segment_generic(za, zb, u,
            (out,uvec,z,α,γ,δ,ε,q)->f_u_inplace_f64!(out,uvec,z,α,γ,δ,ε,q),
            (J,uvec,z,α,γ,δ,ε,q)->jac_u_f64!(J,uvec,z,α,γ,δ,ε,q),
            α, γ, δ, ε, q;
            reltol=reltol, abstol=abstol)
        # renormalize if requested (remove scale from first component)
        if renorm
            scale = u_end[1]
            if !isfinite(real(scale)) || !isfinite(imag(scale)) || scale == 0.0+0.0im
                throw(ErrorException("Non-finite or zero renormalization scale at segment $k"))
            end
            M *= scale
            u = u_end ./ scale
        else
            u = u_end
        end
    end

    return u, M, λ, nodes[end]  # u normalized, multiplier, lambda, final z
end

# ---------- integrate y-system from z0 to zT (starts from given y_state) ----------
function integrate_y_from_initial(y0_vec::Vector{ComplexF64}, z0::ComplexF64, zT::ComplexF64,
                                  γ, δ, ε, α, q; nseg::Int=200, imag_eps=0.0,
                                  reltol=1e-12, abstol=1e-14)
    # path nodes (straight line, possible imag_eps)
    zend = ComplexF64(real(zT), imag_eps)
    nodes = [ z0 + (k/nseg)*(zend - z0) for k in 0:nseg ]

    y = copy(y0_vec)
    for k in 1:nseg
        za = nodes[k]; zb = nodes[k+1]
        y_end = integrate_segment_generic(za, zb, y,
            (out,uvec,z,γ,δ,ε,α,q)->f_y_inplace_f64!(out,uvec,z,γ,δ,ε,α,q),
            (J,uvec,z,γ,δ,ε,α,q)->jac_y_f64!(J,uvec,z,γ,δ,ε,α,q),
            γ, δ, ε, α, q; reltol=reltol, abstol=abstol)
        y = y_end
    end
    return y
end

# ---------- Top-level function: integrate u to z_switch then switch to y for further points ----------
"""
solve_heun_switch(q, α, γ, δ, ε, z_points; z_switch=-10.0, nseg_u=300, nseg_y=200,
                  imag_eps=1e-13, renorm=true, z0_abs=1e-9, reltol=1e-12, abstol=1e-14)

Integrate transformed u from near 0 to z_switch, then switch to original y-equation for z < z_switch.
Works for ComplexF64 inputs and real or complex z_points. Returns Vector{ComplexF64} y-values at z_points.
"""
function solve_heun_switch(q, α, γ, δ, ε, z_points; z_switch=-10.0,
                           nseg_u::Int=300, nseg_y::Int=200,
                           imag_eps::Float64=1e-13, renorm::Bool=true,
                           z0_abs::Float64=1e-9, reltol=1e-12, abstol=1e-14)

    # ensure complex params
    qf = ComplexF64(q); αf = ComplexF64(α); γf = ComplexF64(γ)
    δf = ComplexF64(δ); εf = ComplexF64(ε)

    # prepare output
    results = ComplexF64[]
    # process each z_point independently (safe, simple). Could be optimized.
    for zT_raw in z_points
        zT = ComplexF64(real(zT_raw), imag(zT_raw))
        # if target is between 0 and z_switch (i.e. abs(zT) <= abs(z_switch)), integrate u to zT and return y
        if real(zT) >= z_switch    # assumes z_switch negative; this picks less-negative targets
            # integrate u directly to zT
            try
                u_norm, M, λ, z_final = integrate_u_with_renorm(qf, αf, γf, δf, εf, zT;
                                                              nseg=nseg_u, imag_eps=imag_eps,
                                                              z0_abs=z0_abs, reltol=reltol, abstol=abstol, renorm=renorm)
                # reconstruct y = M * exp(λ z) * u_norm[1]
                y_val = M * exp(λ * ComplexF64(real(zT), imag_eps)) * u_norm[1]
            catch err
                @warn "u-phase integration failed for z=$zT: $err"
                y_val = 0.0 + 0.0im
            end
            # sanitize
            if !isfinite(real(y_val)) || !isfinite(imag(y_val))
                y_val = 0.0 + 0.0im
            end
            push!(results, y_val)
        else
            # target beyond z_switch (more negative): integrate u to z_switch then switch to y-form
            try
                z_switch_c = ComplexF64(z_switch, imag_eps)
                u_norm, M, λ, _ = integrate_u_with_renorm(qf, αf, γf, δf, εf, z_switch_c;
                                                         nseg=nseg_u, imag_eps=imag_eps, z0_abs=z0_abs,
                                                         reltol=reltol, abstol=abstol, renorm=renorm)
                # reconstruct y and y' at z_switch
                y_switch = M * exp(λ * z_switch_c) * u_norm[1]
                yprime_switch = M * exp(λ * z_switch_c) * (λ * u_norm[1] + u_norm[2])
                y0_vec = ComplexF64[y_switch, yprime_switch]
                # integrate y-system from z_switch to zT
                y_end = integrate_y_from_initial(y0_vec, z_switch_c, zT, γf, δf, εf, αf, qf;
                                                 nseg=nseg_y, imag_eps=imag_eps, reltol=reltol, abstol=abstol)
                y_val = y_end[1]
            catch err
                @warn "switch-phase integration failed for target $zT: $err"
                y_val = 0.0 + 0.0im
            end
            if !isfinite(real(y_val)) || !isfinite(imag(y_val))
                y_val = 0.0 + 0.0im
            end
            push!(results, y_val)
        end
    end

    return results
end

# ---------- Example usage ----------
if abspath(PROGRAM_FILE) == @__FILE__
    q = 1.997 + 0.059im
    α = -0.002 + 4.7e-5im
    γ = 1.00 + 1.86im
    δ = 0.999 - 1.984im
    ε = -0.0007 + 0.0im

    z_points = [-1e-6, -0.001, -0.005, -0.01, -1.0, -10.0, -50.0, -200.0, -2000.0]
    ys = solve_heun_switch(q, α, γ, δ, ε, z_points; z_switch=-10.0,
                           nseg_u=300, nseg_y=400, imag_eps=1e-13, renorm=true,
                           z0_abs=1e-9, reltol=1e-12, abstol=1e-14)

    for (z,y) in zip(z_points, ys)
        println("z = $z  y = $y")
    end
end
