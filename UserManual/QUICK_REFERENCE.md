# AxionSR Quick Reference Card

## Installation (2 minutes)

```bash
git clone https://github.com/your-repo/AxionSR.jl.git
cd AxionSR.jl
julia --project
julia> Pkg.instantiate()
```

## Three Essential Functions

### 1. Quick Simulation
```julia
using AxionSR

a_final, M_final, info = super_rad_check(
    M_BH = 5.0,        # Black hole mass [M_⊙]
    aBH = 0.8,         # Initial spin
    massB = 1e-19,     # Axion mass [eV]
    f_a = 1e16,        # Decay constant [GeV]
    age = 1e9          # Age [years]
)
```

### 2. Eigenfrequency
```julia
omega = solve_radial(0.1, M_BH=5.0, a=0.8, n=2, l=1, m=1)
# ω_R = real(omega), ω_I = imag(omega)
```

### 3. Growth Rates
```julia
rates, _, dict = compute_sr_rates(
    [(2,1,1), (3,2,1)],
    M_BH=5.0, aBH=0.8, alpha=0.1
)
```

## Key Parameters

| Parameter | Symbol | Range | Example |
|-----------|--------|-------|---------|
| Black hole mass | $M_{\text{BH}}$ | 0.1–10¹⁰ M☉ | 5.0 |
| Black hole spin | $a$ | 0–0.998 | 0.95 |
| Axion mass | $m_a$ | 10⁻²¹–10⁻¹⁰ eV | 10⁻¹⁹ |
| Decay constant | $f_a$ | 10¹⁴–10¹⁸ GeV | 10¹⁶ |
| Age | $t$ | 10⁶–10¹¹ years | 10⁹ |
| Coupling | $\alpha = Gm_a M$ | 0.001–100 | 0.1 |

## Physics at a Glance

**Superradiance Condition**: $0 < \omega < m \Omega_H$
- Growth rate: $\Gamma = 2|\omega_I|$ [yr⁻¹]
- Timescale: $\tau \sim 1/\Gamma$ (typically 10⁴–10⁸ years)

**Saturation**: Cloud binds energy $\sim M c^2$; bosenova decays cloud

**Spin-down**: $\Delta a \approx E_{\text{cloud}} / M^2$

## Common Workflows

### Test if superradiance affects a black hole
```julia
a_pred, _, _ = super_rad_check(M, a_obs, m_a, f_a, age)
if abs(a_pred - a_obs) < σ_a
    println("Consistent with observation")
else
    println("Axion mass ruled out")
end
```

### Scan parameter space
```julia
for m_a in 10.0 .^ (-20:0.5:-10)
    for f_a in 10.0 .^ (14:0.5:18)
        a_pred, _, _ = super_rad_check(M, a_init, m_a, f_a, age)
        if is_consistent(a_pred, a_obs, σ_a)
            println("$(m_a), $(f_a) allowed")
        end
    end
end
```

### Extract evolution history
```julia
a_final, M_final, info = super_rad_check(...)
t_peak = info["times"][argmax(info["energies"])]
println("Peak cloud energy at t = $(t_peak/1e6) Myr")
```

## Plotting

```julia
using Plots
plot(info["times"]/1e6, info["spins"],
     xlabel="Time [Myr]", ylabel="Spin a", legend=false)
savefig("spin_evolution.pdf")
```

## Key Numbers to Remember

| Quantity | Value | Notes |
|----------|-------|-------|
| Growth rate for $\alpha=0.1$ | 10⁻⁴ yr⁻¹ | Typical superradiant |
| Saturation time (stellar-mass) | 10⁴–10⁶ yr | Rapid on astrophysical timescales |
| Cloud size | $\lambda_C = 1/(m_a c)$ | 10⁻⁵ cm for $m_a = 10⁻²⁰$ eV |
| Bosenova GW luminosity | $\sim 10^{50}$ erg | Observable by LIGO for nearby sources |
| Superradiance threshold | $\alpha \sim 0.1$ | When $\omega_I$ becomes significant |

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Solver doesn't converge | Reduce `Nmax` or increase `K_max` |
| Evolution times out | Use shorter `tau_max` or fewer states |
| NaN/Inf results | Check input parameters are physical |
| Slow performance | Use `interpolate=true` (default) |
| Can't find rate files | Run `Pkg.instantiate()` |

## File Locations

```
AxionSR.jl/
├── src/              # Source code
├── rate_sve/         # Pre-computed rates (100 MB)
├── BH_data/          # Black hole observations
├── test/             # Test suite (44 tests)
└── UserManual/       # This documentation
```

## Mathematical Notation

| Symbol | Meaning |
|--------|---------|
| $\omega = \omega_R + i\omega_I$ | Eigenfrequency (complex) |
| $\Gamma = 2\|\omega_I\|$ | Growth rate [yr⁻¹] |
| $(n,l,m)$ | Quantum numbers (principal, orbital, azimuthal) |
| $\alpha = GM_a M_{\text{BH}}$ | Dimensionless coupling |
| $E_{nlm}(t)$ | Energy in quantum state |
| $a(t)$ | Black hole spin (dimensionless) |

## Important Limits

**Small coupling** ($\alpha < 0.01$):
- Superradiance negligible
- Final spin $\approx$ initial spin

**Strong coupling** ($\alpha > 1$):
- Rapid growth over $\sim 10^5$ years
- Significant spin-down expected
- Nonlinear effects important

**Critical coupling** ($\alpha \sim 0.1$):
- Transition regime
- Growth rates largest
- Most sensitive to observations

## Citation

```bibtex
@manual{AxionSR2025,
  author = {Witte, Samuel D.},
  title = {AxionSR: Computational Framework for Axion Superradiance},
  year = {2025}
}
```

## Getting Help

- **Documentation**: See `UserManual/main.pdf` (full manual)
- **Troubleshooting**: Appendix C of manual
- **Examples**: Section 4 of manual
- **Issues**: GitHub repository issue tracker
- **Tests**: Run `julia test/runtests.jl` to verify installation

## Key Equations

**Superradiance threshold**: $a > a_{\text{threshold}}(n,l,m,\alpha)$ where $\omega_I = 0$

**Growth timescale**: $\tau = 1/(2|\omega_I|) \sim 10^{4-8}$ years

**Spin-down rate**: $\dot{a} \approx -\frac{E_{\text{cloud}}}{M^2}$

**Bosenova condition**: $E_{\text{cloud}} > E_{\text{bind}} \Rightarrow$ condensate decays

---

**Print this card** for quick reference while reading code or papers!

**Last Updated**: November 2025
**Manual Version**: 1.0
