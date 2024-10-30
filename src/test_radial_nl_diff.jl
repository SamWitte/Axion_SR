include("solve_sr_rates.jl")

alp = 0.71586206

M=22.2
a=0.95
mu = alp ./ (GNew .* M)

n1=5
l1=2
m1=2

n2=5
l2=2
m2=2

n3=3
l3=2
m3=2

to_inf = true

rpts=500000
Npts_Bnd=50000
Ntot_safe=5000
eps_fac= BigFloat(1e-3)
Nang=300000
NON_REL = false ## quick approx radial WF and energies
h_mve=10

m=0
l=0

print("Params for Solving non-lin KG: \n")
print("Alpha \t", alp, "\n")
print("spin \t", a, "\n")
print("rpts \t", rpts, "\n")
print("eps_fac \t", eps_fac, "\n")
print("NON_REL \t", NON_REL, "\n")
print("l \t", l, "\n")
print("H mve \t", h_mve, "\n")
print("Rate \t", n1, l1, m1, "\t", n2,l2,m2, "\t", n3,l3,m3, "\n")


@time gf_radial(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts, Npts_Bnd=Npts_Bnd, debug=true,  eps_fac=eps_fac, Ntot_safe=Ntot_safe, m=m, l=l, Nang=Nang, NON_REL=NON_REL,  h_mve=h_mve, to_inf=to_inf)

