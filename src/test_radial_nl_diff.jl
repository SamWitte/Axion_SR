include("solve_sr_rates.jl")

aa = 0.1
# mu = 3e-13

M=22.2
a=0.9
mu = aa ./ (GNew .* M)

n1=3
l1=2
m1=2

n2=3
l2=2
m2=2

n3=2
l3=1
m3=1

to_inf = true

rpts=500000
Npts_Bnd=50000
Ntot_safe=5000
eps_fac= BigFloat(1e-3)
Nang=300000
NON_REL = true ## quick approx radial WF and energies
h_mve=10

m=0
l=0

print("Params for Solving non-lin KG: \n")
print("Alpha \t", GNew .* M .* mu, "\n")
print("spin \t", a, "\n")
print("rpts \t", rpts, "\n")
print("eps_fac \t", eps_fac, "\n")
print("NON_REL \t", NON_REL, "\n")
print("l \t", l, "\n")
print("H mve \t", h_mve, "\n")
print("Rate \t", n1, l1, m1, "\t", n2,l2,m2, "\t", n3,l3,m3, "\n")

# @time integrate_radialEq_2(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts, Npts_Bnd=Npts_Bnd, debug=true,  eps_fac=eps_fac, Ntot_safe=Ntot_safe, m=m, l=l, Nang=Nang, NON_REL=NON_REL,  h_mve=h_mve)

# @time integrate_radialEq_3(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts, Npts_Bnd=Npts_Bnd, , debug=true, sve_for_test=sve_for_test, add_source=add_source, eps_fac=eps_fac, Ntot_safe=Ntot_safe, m=m, l=l, Nang=Nang, NON_REL=NON_REL,  h_mve=h_mve)

@time gf_radial(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts, Npts_Bnd=Npts_Bnd, debug=true,  eps_fac=eps_fac, Ntot_safe=Ntot_safe, m=m, l=l, Nang=Nang, NON_REL=NON_REL,  h_mve=h_mve, to_inf=to_inf)

