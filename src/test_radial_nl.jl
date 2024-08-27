include("solve_sr_rates.jl")
mu = 1e-12
rmaxT = 20

out = direct_solve_radialL1(mu, 22.2, 0.9, 2, 1, 1, 2, 1, 1, 3, 2, 2, debug=true, rpts=1000, Npts_Bnd=4000, tag="_1000_")
print(out, "\n")

out = direct_solve_radialL1(mu, 22.2, 0.9, 2, 1, 1, 2, 1, 1, 3, 2, 2, debug=true, rpts=2000, Npts_Bnd=4000, tag="_1000_")
print(out, "\n")

out = direct_solve_radialL1(mu, 22.2, 0.9, 2, 1, 1, 2, 1, 1, 3, 2, 2, debug=true, rpts=3000, Npts_Bnd=4000, tag="_1000_")
print(out, "\n")

out = direct_solve_radialL1(mu, 22.2, 0.9, 2, 1, 1, 2, 1, 1, 3, 2, 2, debug=true, rpts=4000, Npts_Bnd=4000, tag="_1000_")
print(out, "\n")
