include("solve_sr_rates.jl")

direct_solve_radialL1(mu , 22.2, 0.9, 2, 1, 1, 2, 1, 1, 3, 2, 2, debug=true, rpts=1000, Npts_Bnd=4000)

direct_solve_radialL1(mu , 22.2, 0.9, 2, 1, 1, 2, 1, 1, 3, 2, 2, debug=true, rpts=2000, Npts_Bnd=4000)

direct_solve_radialL1(mu , 22.2, 0.9, 2, 1, 1, 2, 1, 1, 3, 2, 2, debug=true, rpts=3000, Npts_Bnd=4000)

direct_solve_radialL1(mu , 22.2, 0.9, 2, 1, 1, 2, 1, 1, 3, 2, 2, debug=true, rpts=4000, Npts_Bnd=4000)