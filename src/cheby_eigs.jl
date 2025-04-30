using LinearAlgebra
using SpecialFunctions
using WignerSymbols
using HypergeometricFunctions
include("Constants.jl")
# include("heunc.jl")
include("solve_sr_rates.jl")


function eigensys_Cheby(M, atilde, mu, n, l0, m; prec=100, L=4, Npoints=60, Iter=10, debug=false, return_wf=false, der_acc=1e-6, cvg_acc=1e-4, Npts_r=200, nu_guess=nothing, return_nu=false)
    # L field spherical harmonic truncation l-eigenstate
    # Npoints number of Chebyshev interpolation points
    # Iter number of iterations for the non-linear inversion
    
    setprecision(BigFloat, prec)  # Higher than requested precision for calculations
    
    a = BigFloat(atilde) * M  # BH spin
    alph = BigFloat(mu .* GNew)
    M = BigFloat(M)
    n0 = n - l0 - 1 # field overtone number
   
    
    # Define auxiliary quantities
    rminus = M - sqrt(M^2 - a^2)
    rplus = M + sqrt(M^2 - a^2)

    # Helper function for Kronecker delta
    function kronecker_delta(i, j)
        i == j ? 1 : 0
    end

    # Defining ζ-r map
    function ζmap(r)
        return (r - sqrt(4*rplus*(r - rminus) + rminus^2))/(r - rminus)
    end

    function rmap(ζ)
        return (-rminus + 4*rplus + rminus*ζ^2)/(-1 + ζ)^2
    end

    # Black hole function
    function Δ(r)
        return r^2 - 2*M*r + a^2
    end


    # Define coupling function
    function c_coupling(l, j)
        if j-2 <= l <= j+2 && -l <= m <= l && -j <= m <= j
            return (1/3) * kronecker_delta(l, j) +
                   (2/3) * sqrt((2j + 1)/(2l + 1)) *
                   clebschgordan(j, m, 2, 0, l, m) *
                   clebschgordan(j, 0, 2, 0, l, 0)
        else
            return (1/3) * kronecker_delta(l, j)
        end
    end

    # Derivatives of the ζ-r map
    function Dζmap(r)
        return (rminus^2 + 2*r*rplus -
               rminus*(2*rplus + sqrt(rminus^2 + 4*r*rplus - 4*rminus*rplus)))/
               ((r - rminus)^2 * sqrt(rminus^2 + 4*r*rplus - 4*rminus*rplus))
    end

    function DDζmap(r)
        return 1/(r - rminus)^3 * 2 * (r +
               (2*(r - rminus)^2 * rplus^2)/(rminus^2 + 4*r*rplus - 4*rminus*rplus)^(3/2) -
               sqrt(rminus^2 + 4*r*rplus - 4*rminus*rplus) -
               (r - rminus) * (1 - (2*rplus)/sqrt(rminus^2 + 4*r*rplus - 4*rminus*rplus)))
    end

    # Define Chebyshev interpolation points
    ζ = [cos((BigFloat(π) * (2*k + 1))/(2 * (Npoints + 1))) for k in 0:Npoints]  # Chebyshev nodes
    w = [sin((2*BigFloat(π)*k*(Npoints + 2) + π)/(2 * (Npoints + 1))) for k in 0:Npoints]  # Chebyshev weights

    # Derivative matrices in polynomial base
    DerivP = zeros(BigFloat, Npoints+1, Npoints+1)
    for k in 0:Npoints, n in 0:Npoints
        if k != n
            DerivP[k+1, n+1] = w[k+1]/w[n+1]/(ζ[n+1] - ζ[k+1])
        else
            sum_val = sum(w[kk+1]/w[n+1]/(ζ[n+1] - ζ[kk+1]) for kk in 0:Npoints if kk != n)
            DerivP[k+1, n+1] = -sum_val
        end
    end

    # Second derivative matrix
    DDerivP = zeros(BigFloat, Npoints+1, Npoints+1)
    for k in 0:Npoints, n in 0:Npoints
        if k != n
            DDerivP[k+1, n+1] = 2 * DerivP[k+1, n+1] * DerivP[n+1, n+1] -
                               (2 * DerivP[k+1, n+1])/(ζ[n+1] - ζ[k+1])
        else
            sum_val = sum(2 * DerivP[kk+1, n+1]/(ζ[n+1] - ζ[kk+1]) for kk in 0:Npoints if kk != n)
            DDerivP[k+1, n+1] = 2 * DerivP[k+1, n+1] * DerivP[n+1, n+1] + sum_val
        end
    end

    function generalized_laguerre(n, α, x)
        return (-1).^ n ./ factorial(big(n)) .* HypergeometricFunctions.U.(-n, α+1 ,x)
        # return binomial(n + α, n) .* HypergeometricFunctions.M.(-n, α+1 , x)
    end


    # Function for hydrogenic frequency parameter
    function calc_Nν_initial()
        return l0 + 1 + n0 +
            im * (alph * M)^(4*l0 + 2) * ((a*m)/M - 2*alph*rplus) *
            (2^(4*l0 + 2) * factorial(big(2*l0 + 1 + n0)))/
            ((l0 + 1 + n0)^(2*l0 + 1) * factorial(big(n0))) *
            (factorial(big(l0))/(factorial(big(2*l0)) * factorial(big(2*l0 + 1))))^2 *
            prod(j^2 * (1 - a^2/M^2) + ((a*m)/M - 2*alph*rplus)^2 for j in 1:l0)
    end
    
    function calc_Nν_initial_2()
        wI = sr_rates(n, l0, m, mu, M, a)
        wR = GNew .* M .* ergL(n, l0, m, mu, M, a; full=false)
        nuR = (alph .* M).^2 ./ sqrt.((alph .* M).^2 .- wR.^2)
        nuI = nuR.^3 ./ (alph .* M)^3 .* (wI .* GNew .* M) .* sqrt.((nuR.^2 .- (alph .* M)^2) ./ nuR.^2)
        return nuR .+ im * nuI
    end
    

    # Hydrogenic frequency parameter ν (initial value)
    if isnothing(nu_guess)
        Nν = calc_Nν_initial_2()
    else
        Nν = BigFloat(real(nu_guess)) .+ BigFloat(imag(nu_guess))
    end

    # Now we can define ω
    function calc_ω(ν_val)
        return alph * sqrt(1 - (alph^2*M^2)/ν_val^2)
    end
    
    if debug
        println("starting erg \t", calc_ω(Nν).* M)
    end

    # Define P and A values as functions of ν
    function calc_Pplus(ν_val)
        ω_val = calc_ω(ν_val)
        return (m*a - 2*ω_val*M*rplus)/(rplus - rminus)
    end

    function calc_Pminus(ν_val)
        ω_val = calc_ω(ν_val)
        return (m*a - 2*ω_val*M*rminus)/(rplus - rminus)
    end

    function calc_Aplus(ν_val)
        ω_val = calc_ω(ν_val)
        return 2*((m*a - 2*ω_val*M^2)/(rplus - rminus))^2 -
               4*ω_val^2*M^2 + (alph^2 - 2*ω_val^2)*M*(rplus - rminus) +
               (2*M^2 - a^2)*(alph^2 - ω_val^2)
    end

    function calc_Aminus(ν_val)
        ω_val = calc_ω(ν_val)
        return 2*((m*a - 2*ω_val*M^2)/(rplus - rminus))^2 -
               4*ω_val^2*M^2 - (alph^2 - 2*ω_val^2)*M*(rplus - rminus) +
               (2*M^2 - a^2)*(alph^2 - ω_val^2)
    end

    # Define asymptotic behavior functions as functions of ν
    function F(r, ν_val)
        ω_val = calc_ω(ν_val)
        Pplus_val = calc_Pplus(ν_val)
        return ((r - rplus)/(r - rminus))^(im * Pplus_val) *
               (r - rminus)^(-1 + ν_val - (2*alph^2*M^2)/ν_val) *
               exp(-(alph^2*M)/ν_val * (r - rplus))
    end

    function DF(r, ν_val)
        ω_val = calc_ω(ν_val)
        Pplus_val = calc_Pplus(ν_val)
        return -(1/ν_val) * exp(-((M*(r - rplus)*alph^2)/ν_val)) *
               (r - rminus)^(-3 - (2*M^2*alph^2)/ν_val + ν_val) *
               ((r - rplus)/(r - rminus))^(-1 + im*Pplus_val) *
               (2*M^2*(r - rplus)*alph^2 + M*(r - rminus)*(r - rplus)*alph^2 +
                ν_val*(r + im*Pplus_val*(rminus - rplus) + rplus*(-1 + ν_val) - r*ν_val))
    end

    function DDF(r, ν_val)
        ω_val = calc_ω(ν_val)
        Pplus_val = calc_Pplus(ν_val)
        return (1/((r - rplus)^2 * ν_val^2)) * exp(-((M*(r - rplus)*alph^2)/ν_val)) *
               (r - rminus)^(-3 - (2*M^2*alph^2)/ν_val + ν_val) *
               ((r - rplus)/(r - rminus))^(im*Pplus_val) *
               (4*M^4*(r - rplus)^2*alph^4 + 4*M^3*(r - rminus)*(r - rplus)^2*alph^4 -
                2*M*(r - rminus)*(r - rplus)*alph^2*ν_val*(-im*Pplus_val*(rminus - rplus) +
                rplus + r*(-1 + ν_val) - rplus*ν_val) + M^2*(r - rplus)*alph^2*(r^3*alph^2 -
                rminus^2*rplus*alph^2 - r^2*(2*rminus + rplus)*alph^2 + 4*im*Pplus_val*rminus*ν_val +
                2*rplus*ν_val*(-3 - 2*im*Pplus_val + 2*ν_val) + r*(rminus^2*alph^2 + 2*rminus*rplus*alph^2 +
                2*(3 - 2*ν_val)*ν_val)) + ν_val^2*(-Pplus_val^2*(rminus - rplus)^2 -
                im*Pplus_val*(rminus - rplus)*(rminus + 3*rplus - 2*rplus*ν_val) +
                2*r*(-2 + ν_val)*(-im*Pplus_val*(rminus - rplus) + rplus - rplus*ν_val) +
                r^2*(2 - 3*ν_val + ν_val^2) + rplus^2*(2 - 3*ν_val + ν_val^2)))
    end

    # Hydrogenic radial wave function
    function R(r, ν_val)
        return ((2*r*M*alph^2)/(l0 + n0 + 1))^l0 *
               exp(-((r*M*alph^2)/(l0 + n0 + 1))) *
               generalized_laguerre(n - l0 - 1, 2 * l0 + 1, (2*r*M*alph^2)/(l0 + n0 + 1))
               
    end

    # Function to build equation coefficients
    function build_coefficients(ν_val)
        ω_val = calc_ω(ν_val)
        Pplus_val = calc_Pplus(ν_val)
        Pminus_val = calc_Pminus(ν_val)
        Aplus_val = calc_Aplus(ν_val)
        Aminus_val = calc_Aminus(ν_val)
        
        C1 = zeros(Complex{BigFloat}, L+1, Npoints+1)
        C2 = zeros(Complex{BigFloat}, L+1, Npoints+1)
        C3plus = zeros(Complex{BigFloat}, L+1, Npoints+1)
        C3minus = zeros(Complex{BigFloat}, L+1, Npoints+1)
        
        for l in 0:L, n in 0:Npoints
            r_n = rmap(ζ[n+1])
            C1[l+1, n+1] = (1/(r_n - rplus) + 1/(r_n - rminus)) * 1/Dζmap(r_n) +
                           (2*DF(r_n, ν_val))/F(r_n, ν_val) * 1/Dζmap(r_n) +
                           DDζmap(r_n)/(Dζmap(r_n)^2)
            
            C2[l+1, n+1] = 1/(Dζmap(r_n))^2 * DDF(r_n, ν_val)/F(r_n, ν_val) +
                           (1/(r_n - rplus) + 1/(r_n - rminus)) *
                           1/(Dζmap(r_n))^2 * DF(r_n, ν_val)/F(r_n, ν_val) +
                           1/(Dζmap(r_n))^2 * (Pplus_val^2/(r_n - rplus)^2 +
                           Pminus_val^2/(r_n - rminus)^2 -
                           Aplus_val/((rplus - rminus)*(r_n - rplus)) +
                           Aminus_val/((rplus - rminus)*(r_n - rminus)) - (alph^2 - ω_val^2)) -
                           1/(Dζmap(r_n))^2 * 1/Δ(r_n) *
                           (l*(l + 1) + a^2*c_coupling(l, l)*(alph^2 - ω_val^2))
            
            C3plus[l+1, n+1] = -(1/(Dζmap(r_n))^2) * 1/Δ(r_n) *
                               a^2 * (alph^2 - ω_val^2) * c_coupling(l, l + 2)
            
            C3minus[l+1, n+1] = -(1/(Dζmap(r_n))^2) * 1/Δ(r_n) *
                                a^2 * (alph^2 - ω_val^2) * c_coupling(l, l - 2)
        end
        
        return C1, C2, C3plus, C3minus
    end

    # Function to build the system matrix
    function build_matrix(ν_val)
        C1, C2, C3plus, C3minus = build_coefficients(ν_val)
        
        # Size of the full matrix
        N = (L+1)*(Npoints+1)
        MatrixOp = zeros(Complex{BigFloat}, N, N)
        
        for l in 0:L
            for n in 0:Npoints
                row = l*(Npoints+1) + n + 1
                
                # Diagonal part (l,l)
                for k in 0:Npoints
                    col = l*(Npoints+1) + k + 1
                    MatrixOp[row, col] = DDerivP[k+1, n+1] + C1[l+1, n+1]*DerivP[k+1, n+1]
                    if k == n
                        MatrixOp[row, col] += C2[l+1, n+1]
                    end
                end
                
                # Off-diagonal part (l,l-2)
                if l >= 2
                    col = (l-2)*(Npoints+1) + n + 1
                    MatrixOp[row, col] = C3minus[l+1, n+1]
                end
                
                # Off-diagonal part (l,l+2)
                if l+2 <= L
                    col = (l+2)*(Npoints+1) + n + 1
                    MatrixOp[row, col] = C3plus[l+1, n+1]
                end
            end
        end
        
        return MatrixOp
    end

    # Function to build the derivative of the system matrix with respect to ν
    function build_derivative_matrix(ν_val)
        # This is a complex operation, would need symbolic differentiation
        # Placeholder implementation - actual implementation would need careful derivative calculations
        
        # Use numerical differentiation as an approximation
        ε = der_acc
        return (build_matrix(ν_val .* (1 .+ ε)) - build_matrix(ν_val .* (1 .- ε)))/(2*ν_val*ε)
    end

    # Initialize solution vectors for non-linear inverse iteration
    Nν_values = Vector{Complex{BigFloat}}(undef, Iter+1)
    Nν_values[1] = Nν  # Initial value

    # Function to build initial guess vector
    function build_initial_vector()
        N = (L+1)*(Npoints+1)
        v = zeros(Complex{BigFloat}, N)
        
        counter = 1
        for j in 0:L
            for k in 0:Npoints
                v[counter] = (kronecker_delta(j, l0) * R(rmap(ζ[k+1]), Nν))/F(rmap(ζ[k+1]), Nν)
                counter += 1
            end
        end
        
        return v
    end

    # Normalization vector for iteration
    function build_normalization_vector()
        N = (L+1)*(Npoints+1)
        return [exp(im * jk * π/N)/sqrt(N) for jk in 1:N]
    end

    # Perform non-linear inverse iteration
    function run_nonlinear_inverse_iteration()
        # Initialize vectors
        Bnum = Vector{Vector{Complex{BigFloat}}}(undef, Iter+1)
        BnumNorm = Vector{Vector{Complex{BigFloat}}}(undef, Iter+1)
        
        # Initial guess
        Bnum[1] = build_initial_vector()
        BnumNorm[1] = Bnum[1] / sqrt(dot(conj(Bnum[1]), Bnum[1]))
        
        # Normalization vector
        VEC = build_normalization_vector()
        
        # Iteration
        final_idx = 1
        i = 1
        while i < (Iter + 1)
            A_i = build_matrix(Nν_values[i])
            b_i = build_derivative_matrix(Nν_values[i]) * BnumNorm[i]
            
            # Solve linear system
            Bnum[i+1] = A_i \ b_i
            
            # Update eigenvalue estimate
            Nν_values[i+1] = Nν_values[i] - dot(conj(VEC), BnumNorm[i]) / dot(conj(VEC), Bnum[i+1])
            
            # Normalize solution vector
            BnumNorm[i+1] = Bnum[i+1] / sqrt(dot(conj(Bnum[i+1]), Bnum[i+1]))
            
            if debug
                println("Iteration $i: ν = $(Nν_values[i+1])")
            end
            
            if i > 1
                erg_test_init = calc_ω(Nν_values[i])
                erg_test = calc_ω(Nν_values[i+1])
                test_convR = abs.(real.(erg_test .- erg_test_init) ./ real.(erg_test_init))
                test_convI = abs.(imag.(erg_test .- erg_test_init) ./ imag.(erg_test_init))
                if (test_convR < cvg_acc)&&(test_convI < cvg_acc)
                    final_idx = i
                    i = Iter + 1
                end
            end
            i += 1
            
            
        end
        
        return Nν_values, Bnum, BnumNorm, final_idx
    end


    # Run the algorithm
    Nν_values, Bnum, BnumNorm, final_idx = run_nonlinear_inverse_iteration()
    erg_out = calc_ω(Nν_values[final_idx]) .* M
    

    
    if debug
        println("Final ν value: ", erg_out)
    end
    
    
    if !return_wf
        if !return_nu
            return real(erg_out), imag(erg_out) # real(G * M * omega), imag(G * M * omega)
        else
            return real(erg_out), imag(erg_out), Nν_values[final_idx] # real(G * M * omega), imag(G * M * omega),
        end
    end

    
    
    x_values = Float64[]
    y_values = Complex[]

##### this it original way
    for j in 1:Npoints
        # Calculate the current r value
        r_val = rmap(ζ[j])

        # Calculate ω based on the final Nν value from iteration
        ω = alph * sqrt(1 - (alph^2 * M^2)/(Nν_values[final_idx])^2)

        # Calculate Pplus for the current ν value
        Pplus = calc_Pplus(Nν_values[final_idx])

        # Calculate the expression from the Mathematica code
        term1 = ((r_val - rplus)/(r_val - rminus))^(im * Pplus)
        term2 = (r_val - rminus)^(-1 + (alph^2 * M)/sqrt(alph^2 - ω^2) - 2 * M * sqrt(alph^2 - ω^2))
        term3 = exp(-sqrt(alph^2 - ω^2) * (r_val - rplus))
        term4 = BnumNorm[final_idx][j + l0 * (Npoints + 1)]

        result = term1 * term2 * term3 * term4
        
        push!(x_values, r_val)
        push!(y_values, result)
    end
    rlist = reverse(x_values) ./ M
    y_values = reverse(y_values)
    
    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
    # norm result
    nm2 = trapz(y_values .* conj.(y_values) .* rlist.^2, rlist)
    y_values ./= sqrt.(nm2)
    
    if debug
        println(rlist)
        println(y_values)
        println("Norm test \t", trapz(y_values .* conj.(y_values) .* rlist.^2, rlist))
    end
 
    if !return_nu
        return real(erg_out), imag(erg_out), rlist, y_values # everything normalized by GM, radial WF may not resolve at large r
    else
        return real(erg_out), imag(erg_out), rlist, y_values, Nν_values[final_idx]
    end
end


# test run of system
# @time eigensys_Cheby(1, 0.95, 0.01 ./ GNew, 4, 3, 3, debug=true, return_wf=true, L=4, Npoints = 60, Iter = 20,  der_acc=1e-6, cvg_acc=1e-3, prec=100, Npts_r=50)
