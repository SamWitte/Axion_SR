include("solve_sr_rates.jl")
using DelimitedFiles

M = 22.2
a=0.9

alpL = LinRange(0.1, 0.3, 20)

n1=2
l1=1
m1=1

n2=2
l2=1
m2=1

n3=3
l3=2
m3=2

rpts=50000
Npts_Bnd=4000
rmaxT=500
Ntot_safe=5000
sve_for_test=false
add_source=true
eps_fac= BigFloat(1e-4)
dtmax_slv = 1e-1
reltol=1e-10
Nang=600000
n_times=300000


m=0
l=0

ratL = zeros(length(alpL), 2)
for i in 1:length(alpL)
    mu = alpL[i] ./ (GNew .* M)
    NR =  integrate_radialEq_2(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts, Npts_Bnd=Npts_Bnd, rmaxT=rmaxT, debug=true, sve_for_test=sve_for_test, add_source=add_source, eps_fac=eps_fac, Ntot_safe=Ntot_safe, dtmax_slv=dtmax_slv, m=m, l=l, reltol=reltol, Nang=Nang, NON_REL=true, n_times=n_times)
    REL =  integrate_radialEq_2(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts, Npts_Bnd=Npts_Bnd, rmaxT=rmaxT, debug=true, sve_for_test=sve_for_test, add_source=add_source, eps_fac=eps_fac, Ntot_safe=Ntot_safe, dtmax_slv=dtmax_slv, m=m, l=l, reltol=reltol, Nang=Nang, NON_REL=false, n_times=n_times)
    
    ratL[i, :] = [alpL[i] REL ./ NR]
    print("\n", "Out \t ", alpL[i], "\t", REL ./ NR, "\n\n")
end

writedlm("test_store/RatioRates_test.dat", ratL)


