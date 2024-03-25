using SpecialFunctions
using Random
using OrdinaryDiffEq
using Statistics
using Distributions
using ForwardDiff: gradient, derivative, Dual, Partials, hessian
using DelimitedFiles
using NLsolve
using DifferentialEquations
include("Constants.jl")


function ergL(n, l, m, massB, MBH)
    alph = GNew .* MBH .* massB
    return massB .* (1.0 .- alph.^2 ./ (2 .* n.^2) - alph.^4 ./ (8 * n.^4) + alph.^4 ./ (8 * n^3) * (-6 / (2 * l + 1) + 2 / n))
end


function sr_rates(n, l, m, massB, MBH, aBH; impose_low_cut=0.01, solve_322=true)
    if (n==3)&&(l==2)&&(m==1)&&(solve_322==false)
        return 0.0
    end
    

    alph = GNew .* MBH .* massB
    
#    if (alph ./ l < impose_low_cut)&&(MBH < 1e2)
#        # We expect binaries to disrupt.
#        return 0.0
#    end
    
    
    rP = nothing
    if aBH .> 0.998
        rP = 1.0 .+ sqrt.(1 - 0.998 .^2)
        aBH = 0.998
    elseif aBH .< 0.0
        rP = 2.0
        aBH = 0.0
    else
        rP = 1.0 .+ sqrt.(1 - aBH.^2)
    end
    rP *= (GNew .* MBH)
    OmegaH = aBH ./ (2 .* (GNew .* MBH) .* (1 .+ sqrt.(1 .- aBH.^2)))
    Anl = 2 .^(4 .* l .+ 1) .* factorial(Int(l .+ n)) ./ (n.^(2 .* l .+ 4) .* factorial(n .- l .- 1))
    Anl *= (factorial(Int(l)) ./ (factorial(Int(2 .* l)) .* factorial(Int(2 .* l .+ 1)))).^2
    Chilm = 1.0
    # erg = ergL(n, l, m, massB, MBH)
    for k in 1:Int(l)
        Chilm *= (k.^2 .* (1.0 .- aBH.^2) .+ (aBH * m .- 2 .* (rP ./ (GNew .* MBH)) .* alph).^2)
        # Chilm *= (k.^2 .* (1.0 .- aBH.^2) .+ (aBH * m .- 2 .* (rP ./ (GNew .* MBH)) .* alph).^2)
    end
    Gamma_nlm = 2 * massB .* rP .* (m .* OmegaH .- ergL(n, l, m, massB, MBH)) .* alph.^(4 .* l + 4) .* Anl .* Chilm
    

    if Gamma_nlm > 0.0
        return Gamma_nlm
    else
        return 0.0
    end
end



function find_im_part(mu, M, a, n, l, m; debug=false, Ntot=200, iter=50, xtol=1e-20)
    # n not used at moment, in reality l ~ n + l?

    OmegaH = a ./ (2 .* (GNew .* M) .* (1 .+ sqrt.(1 .- a.^2)))
    if (ergL(n, l, m, mu, M) < m .* OmegaH)
        
        
        b = sqrt.(1 - a.^2)
        alph = GNew * M * mu
        # w0 = alph .* (1 - alph.^2 ./ (2 .* (n + l + 1).^2)) * (1 + 1e-20 .* im)
        w0 = ergL(n, l, m, mu, M) .* GNew * M * (1 + 1e-20 .* im)
#        guess = sr_rates(n, l, m, mu, M, a)
        
        function wrapper!(F, x)
            wR = x[1]  # real
            wI = x[2]  # imag
            erg = wR + im * wI
            
            q = - sqrt.(alph.^2 .- erg.^2)
            
            gam = im * a * sqrt.(erg.^2 .- alph.^2)
            LLM = l * (l + 1)
            LLM += (-1 + 2 * l * (l + 1) - 2 * m.^2) * gam.^2 ./ (-3 + 4 * l * (l + 1))
            LLM += ((l - m - 1 * (l - m) * (l + m) * (l + m - 1)) ./ ((-3 + 2 * l) * (2 * l - 1).^2) - (l + 1 - m) * (2 * l - m) * (l + m + 1) * (2 + l + m) ./ ((3 + 2 * l).^2 * (5 + 2 * l))) * gam.^4 ./ (2 * (1 + 2 * l))
            LLM += (4 * ((-1 + 4 * m^2) * (l * (1 + l) * (121 + l * (1 + l) * (213 + 8 * l * (1 + l) * (-37 + 10 * l * (1 + l)))) - 2 * l * (1 + l) * (-137 + 56 * l * (1 + l) * (3 + 2 * l * (1 + l))) * m^2 + (705 + 8 * l * (1 + l) * (125 + 18 * l * (1 + l))) * m^4 - 15 * (1 + 46 * m^2))) * gam^6) / ((-5 + 2 * l) * (-3 + 2 * l) * (5 + 2 * l) * (7 + 2 * l) * (-3 + 4 * l * (1 + l))^5)

   
            c0 = 1.0 .- 2.0 * im * erg - 2 * im ./ b .* (erg .- a .* m ./ 2.0)
            c1 = -4.0 .+ 4 * im * (erg - im * q * (1.0 + b)) + 4 * im / b * (erg - a * m / 2.0) .- 2.0 * (erg.^2 .+ q.^2) ./ q
            c2 = 3.0 - 2 * im * erg - 2.0 * (q.^2 - erg.^2) ./ q - 2.0 * im / b * (erg - a * m ./ 2)
            c3 = 2.0 * im * (erg - im * q).^3 ./ q .+ 2 * (erg .- im * q).^2 .* b + q.^2 * a.^2 .+ 2 * im * q * a * m - LLM - 1 - (erg - im * q).^2 ./ q .+ 2 * q * b + 2 * im / b * ( (erg - im * q).^2 ./ q + 1.0) * (erg .- a * m / 2)
            c4 = (erg .- im * q).^4 ./ q.^2 .+ 2 * im * erg * (erg .- im * q).^2 ./ q .- 2 * im ./ b .* (erg .- im * q).^2 ./ q .* (erg .- a * m ./ 2)
            
            function alphaN(nn)
                return nn.^2 .+ (c0 .+ 1) * nn + c0
            end
            function betaN(nn)
                return -2 * nn.^2 .+ (c1 .+ 2) * nn + c3
            end
            function gammaN(nn)
                return nn.^2 .+ (c2 .- 3) * nn + c4
            end
            
            temp = 1
            for i in Ntot:-1:1
                temp = gammaN(i) ./ (betaN(i) .- alphaN(i) .* temp)
            end
            recur = betaN(0) ./ alphaN(0) .- temp
        
            F[1] = real(recur)
            F[2] = imag(recur)
            
            # print(x, "\t", F, "\n")
        end
        
        
        sol = nlsolve(wrapper!, [real(w0), imag(w0)], autodiff = :forward, xtol=xtol, ftol=1e-20, iterations=iter)
        if debug
            print(sol, "\n\n")
            
            print(sol.zero  , "\n")
            Ff = zeros(2)
            wrapper!(Ff, sol.zero)
            print(Ff, "\n")
        end
        if sol.zero[2] > 0
            return sol.zero[2]  # im(G * M * omega)
        else
            return 0.0
        end
    else
        return 0.0
    end
end
