include("solve_sr_rates.jl")

mu = 2.5e-13

M=22.2
a=0.9

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
eps_fac= BigFloat(1e-2)
dtmax_slv = 1e-1
reltol=1e-10
Nang=300000
n_times=300000
NON_REL = true
h_mve=0.1

m=0
l=0

print("Params for Solving non-lin KG: \n")
print("Alpha \t", GNew .* M .* mu, "\n")
print("spin \t", a, "\n")
print("rpts \t", rpts, "\n")
print("rmaxT \t", rmaxT, "\n")
print("add_source \t", add_source, "\n")
print("eps_fac \t", eps_fac, "\n")
print("dtmax_slv \t", dtmax_slv, "\n")
print("reltol \t", reltol, "\n")
print("NON_REL \t", NON_REL, "\n")
print("l \t", l, "\n")


@time integrate_radialEq_2(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts, Npts_Bnd=Npts_Bnd, rmaxT=rmaxT, debug=true, sve_for_test=sve_for_test, add_source=add_source, eps_fac=eps_fac, Ntot_safe=Ntot_safe, dtmax_slv=dtmax_slv, m=m, l=l, reltol=reltol, Nang=Nang, NON_REL=NON_REL, n_times=n_times, h_mve=h_mve)


