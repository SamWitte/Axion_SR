# Phase 2.2: Detailed Merge Plan for solve_system Functions

## Analysis of Differences

### Common Code (80% overlap)
Both functions share:
1. Basic initialization (alph, default_reltol, reltol_Thres)
2. State vector setup (y0, reltol arrays)
3. Tolerance definitions (e_init, def_spin_tol)
4. ODE callback structure (check_timescale, affect_timescale, check_spin, affect_spin)
5. Time span setup (tspan, saveat)
6. ODE solver instantiation and solving
7. Output processing (extracting spinBH, MassB arrays)

### Key Differences

| Aspect | solve_system | solve_system_spinone |
|--------|--------------|----------------------|
| **Quantum Levels** | Nmax up to 8 | Fixed at 1 (spin-1) |
| **Parameters** | 19 parameters | 4 parameters |
| **Level Setup** | Complex nested loops (lines 59-91) | Simple: idx_lvl=1, m_list=[1] |
| **Rate Computation** | `compute_sr_rates()` with interpolation | `precomputed_spin1()` direct |
| **Mvars Structure** | [mu, fa, Emax2, aBH, M_BH, impose_low_cut] (6 items) | [mu, aBH, M_BH] (3 items) |
| **Emax2 Computation** | Complex emax_211() call | Not needed |
| **RHS_ax! Rates** | Loop through SR_rates array | Single precomputed_spin1 call |
| **Spin Check** | More complex spin checking (line 443) | Simpler check (line 560) |

## Merge Strategy

### Step 1: Create Unified Function Signature
Use parameter `spinone::Bool = false` to distinguish modes.

```julia
function solve_system(mu, fa_or_nothing, aBH, M_BH, t_max;
                     spinone::Bool=false, n_times=10000, debug=false,
                     impose_low_cut=0.01, return_all_info=false,
                     eq_threshold=1e-100, stop_on_a=0, abstol=1e-30,
                     non_rel=true, high_p=true, N_pts_interp=200,
                     N_pts_interpL=200, Nmax=3, cheby=true)
```

**Issue**: `fa` parameter only exists in non-spinone mode.
**Solution**: Make it optional or use a wrapper function.

### Step 2: Extract Shared Initialization
Create helper: `function initialize_solver_state(mu, M_BH, aBH, idx_lvl, non_rel, high_p)`

Returns:
- `default_reltol`
- `reltol_Thres`
- `e_init`
- `y0` (initialized with correct size)
- `reltol` (tolerance array)

### Step 3: Extract Quantum Level Setup
Create conditionals:

```julia
if spinone
    idx_lvl = 1
    m_list = [1]
    bn_list = []  # Not used in spinone
    modes = []
    Mvars = [mu, aBH, M_BH]
    # No Emax2 needed
else
    # Current solve_system logic for level setup (lines 59-91)
    # Compute Emax2
    Mvars = [mu, fa, Emax2, aBH, M_BH, impose_low_cut]
end
```

### Step 4: Extract Rate Computation
Create conditional for SR_rates initialization:

```julia
if spinone
    wR, wI = precomputed_spin1(alph, aBH, M_BH)
    SR_rates = [2 .* wI]
else
    SR_rates, interp_funcs, interp_dict = compute_sr_rates(modes, M_BH, aBH, alph, cheby=cheby)
end
```

### Step 5: Create Unified RHS_ax! with Branching
The RHS function has minor differences in spin checking and rate updates.

```julia
function RHS_ax!(du, u, Mvars, t)
    u_real = exp.(u)

    if spinone
        mu, aBH_i, M_BH_i = Mvars
        # spinone-specific spin checking (simpler)
    else
        mu, fa, Emax2, aBH_i, M_BH_i, impose_low_cut = Mvars
        # standard-mode spin checking (more complex)
    end

    # Rest of shared logic
    # ...

    if spinone
        wR, wI = precomputed_spin1(alph, u_real[spinI], u_real[massI])
        SR_rates = [2 .* wI]
    else
        SR_rates = [func(u_real[spinI]) for func in interp_funcs]
    end
end
```

## Implementation Steps

1. **Create backup**: `git branch phase-2-2-merge`
2. **Implement unified function** with spinone parameter
3. **Move spinone-specific code** into conditional blocks
4. **Extract common initialization** into helpers
5. **Test thoroughly**:
   - Run with spinone=false (should match current solve_system)
   - Run with spinone=true (should match current solve_system_spinone)
6. **Keep old functions as deprecated wrappers** (for backward compatibility)
7. **Commit and cleanup**

## Expected Outcome

```
Before:
- solve_system: 456 lines
- solve_system_spinone: 195 lines
- Total: 651 lines, ~80% duplicate

After:
- solve_system: 350-380 lines (unified with branching)
- Improvement: ~150-180 lines removed
```

## Risk Mitigation

- ✅ Extract shared code in small, testable chunks
- ✅ Maintain identical behavior with spinone parameter
- ✅ Keep old functions as wrappers during transition
- ✅ Test both code paths thoroughly
- ✅ Use git branches for safety

## Success Criteria

- [ ] New unified function handles both spinone and standard modes
- [ ] Numerical results identical to pre-refactor code
- [ ] Code size reduced by ~25%
- [ ] All tests pass
- [ ] Backward compatibility maintained
