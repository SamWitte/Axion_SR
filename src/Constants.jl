using ForwardDiff: Dual

const c_km = 2.99792e5
const c_cm = 2.99792e10
const hbar = 6.582119e-16
const GNew = 7484169213.942707 # 1 / (M_SM * eV)
const M_to_eV = 1.1219176708765e+66
const M_pl = 2.435e18 # GeV
const θR = 0.0
const err_CS = 1e-5
const m_elec = 5.11e5
const G_to_eV2 = 1.95e-2
const e_charge = 0.3
const a_fs = 1.0 ./ 137.0

# Seed 3-dim vector with dual partials for gradient calculation
global seed = x -> [map(y -> Dual(y, (1., 0., 0.)), x[:,1]) map(y -> Dual(y, (0., 1., 0.)), x[:,2]) map(y -> Dual(y, (0., 0., 1.)), x[:,3])]
# Extract gradient from dual
global grad = x -> [map(x -> x.partials[1], x) map(x -> x.partials[2], x) map(x -> x.partials[3], x)]


function File_Name_Out(Mass_a, Ax_g, Mns, Rns, v_NS, B0, P, θm, eta_fill, Ntrajs; file_tag="", indx=nothing, mag_mod="Dipole", B_DQ=0.0, θQ=0.0, reflect=false, dead=false, dead_rmax=30.0)
    
    if indx != nothing
        file_tag *= string(indx);
    end
    
    
    fileN = "results/NSrun_Ma_"*string(round(Mass_a, sigdigits=3))*"_AxionG_"*string(Ax_g)
    
    # neutron star physical params
    fileN *= "_B0_"*string(round(B0,sigdigits=4))*"_P_"*string(round(P, sigdigits=4))*"_ThetaM_"*string(round(θm, sigdigits=3));
    if mag_mod != "Dipole"
        fileN *= "_Bdq_"*string(B_DQ)*"_ThetaQ_"*string(θQ);
    end
    fileN *= "_Mns_"*string(Mns)*"_Rns_"*string(Rns);
    fileN *= "_vx_"*string(round(v_NS[1], sigdigits=3));
    fileN *= "_vy_"*string(round(v_NS[2], sigdigits=3));
    fileN *= "_vz_"*string(round(v_NS[3], sigdigits=3));
    
    fileN *= "_etaFill_"*string(eta_fill);
    if reflect
        fileN *= "_reflect_"
    end
    if dead
        fileN *= "_Dead_Rmax_"*string(dead_rmax);
    end
    
    fileN *= "_Ntrajs_"*string(Ntrajs)*"_"*file_tag*"_.h5";

end
