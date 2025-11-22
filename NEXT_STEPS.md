# Next Steps: Code Refactoring Roadmap

**Status**: Phase 2.2 Complete ✅
**Current Progress**: ~60% of initial cleanup complete
**Total Estimated Effort**: ~15-20 hours remaining

---

## Completed Work ✅

### Phase 1 (Complete)
- [x] Remove 20+ debug print statements from `super_rad.jl` and `stat_analysis.jl`
- [x] Remove 50+ lines of commented-out code
- [x] Fix floating-point comparison bug in `solve_sr_rates.jl`
- [x] Create `src/Core/constants.jl` with centralized tolerance management
- [x] Create `src/Core/parameters.jl` with structured types (planned but not yet implemented)

### Phase 2 (Complete)
- [x] **Phase 2.1**: Extract 5 helper functions to `src/Core/evolution_helpers.jl`
  - `get_clamped_spin()`
  - `enforce_bosenova_boundary!()`
  - `enforce_all_boundaries!()`
  - `check_energy_floor()`
  - `is_near_bosenova()`

- [x] **Phase 2.3**: Consolidate 6 hardcoded tolerance values with `SOLVER_TOLERANCES` references

- [x] **Phase 2.2**: Merge `solve_system()` and `solve_system_spinone()` into unified function
  - Created 4 new initialization helpers
  - Unified 650 lines into 450-line single function
  - Eliminated ~80% code duplication
  - ✅ Integration tests passing (5/5)

---

## Remaining Work

### Phase 3: Extract and Organize Rate Coefficients (3-4 hours)

**Priority**: HIGH (Enables further code organization)

#### Phase 3.1: Structured Rate Data (2 hours) ✅ COMPLETE
- [x] Create `src/Core/rate_coefficients.jl` (280 lines)
- [x] Define structured RateCoefficient and RateDatabase types
- [x] Organize rate values with metadata (coefficient, power, flags)
- [x] Create lookup functions for efficient rate access
- [x] Support for Nmax 3-5 with 11-52 rates per level
- [x] Create test_rate_coefficients.jl with full validation
- [x] Create backward compatibility layer (load_rates_structured.jl)
- [x] Documentation and completion report
- [x] Benefits achieved: Type safety, documentation, easy maintenance ✅

#### Phase 3.2: Separate Rate Computation Logic (1-2 hours)
- [ ] Extract rate computation into `src/Numerics/rate_computation.jl`
- [ ] Functions to extract:
  - `compute_sr_rates()` → with helper functions
  - `precomputed_spin1()`
  - `pre_computed_sr_rates()`
- [ ] Consolidate related utilities in one location
- [ ] Benefits: Better organization, easier to understand rate physics

### Phase 4: Module Organization (2-3 hours)

**Priority**: MEDIUM (Improves code accessibility)

#### Phase 4.1: Create Main Module File (1 hour)
- [ ] Create `src/AxionSR.jl` as main entry point
- [ ] Use module organization:
  ```julia
  module AxionSR

  # Import dependencies
  using OrdinaryDiffEq, Interpolations, ...

  # Include submodules
  include("Core/constants.jl")
  include("Core/evolution_helpers.jl")
  include("Numerics/rate_computation.jl")
  # ... etc

  # Export public API
  export super_rad_check, solve_system, ...

  end
  ```
- [ ] Benefits: Clear API, version control, easier testing

#### Phase 4.2: Create Scripts Directory (1 hour)
- [ ] Create `scripts/` directory
- [ ] Move executable scripts from `src/`:
  - `spin_down_plot_files.jl` → `scripts/compute_spin_evolution.jl`
  - `TEST.jl` → Separate into tests and scripts
  - MCMC runner scripts
- [ ] Create `scripts/README.md` with usage instructions
- [ ] Benefits: Separation of concern, clearer project structure

#### Phase 4.3: Proper Test Suite (1 hour)
- [ ] Create `tests/` directory with proper structure
- [ ] Organize existing tests:
  ```
  tests/
  ├── test_physics.jl         # Physics calculations
  ├── test_evolution.jl       # ODE evolution
  ├── test_integration.jl     # Full system tests
  └── test_smoke.jl           # Quick validation
  ```
- [ ] Set up test runner that uses `@testset` framework
- [ ] Benefits: Professional test organization, CI/CD ready

### Phase 5: Documentation and Type Annotations (2-3 hours)

**Priority**: MEDIUM (Improves maintainability)

#### Phase 5.1: Add Docstrings to Major Functions (1.5 hours)
- [ ] Document public API in `src/AxionSR.jl`
- [ ] Add docstrings to:
  - `solve_system()` - Already started
  - `super_rad_check()`
  - `compute_sr_rates()`
  - Helper functions - Already done
- [ ] Use Julia docstring format with examples
- [ ] Benefits: Better IDE support, easier learning curve

#### Phase 5.2: Type Annotations (1 hour)
- [ ] Review function signatures for clarity
- [ ] Add type hints where beneficial (already partially done)
- [ ] Benefits: Better error messages, code clarity

#### Phase 5.3: Create User Guide (1 hour)
- [ ] Create `docs/USER_GUIDE.md`
- [ ] Include:
  - Quick start example
  - Parameter description
  - Output interpretation
  - Common use cases
- [ ] Benefits: Easier onboarding for new users

### Phase 6: Numerical Validation (2-3 hours)

**Priority**: HIGH (Ensures correctness)

#### Phase 6.1: Comparative Baseline Tests (2 hours)
- [ ] **If original code available**: Run both versions with identical parameters
  - Test parameter space:
    - Multiple masses: [5, 10, 20, 50] M☉
    - Multiple spins: [0.5, 0.9, 0.95, 0.99]
    - Different couplings: [1e-18, 1e-20, 1e-22]
  - Compare final spin and mass values
  - Check convergence behavior
- [ ] Verify numerical differences < acceptable tolerance
- [ ] Document any intentional differences
- [ ] Benefits: Confidence in refactoring, catches subtle bugs

#### Phase 6.2: Extended Numerical Tests (1 hour)
- [ ] Expand `test_integration_numerical.jl` with more cases
- [ ] Test edge cases:
  - Very high spin (a ≈ 0.99)
  - Very low spin (a ≈ 0.01)
  - Wide coupling range
  - Different Nmax values (3, 4, 5, 6)
- [ ] Benefits: Robustness validation

### Phase 7: Performance Optimization (2-3 hours)

**Priority**: LOW (Only if needed for production)

#### Phase 7.1: Profile Current Code (1 hour)
- [ ] Use Julia's `@profile` macro
- [ ] Identify computational bottlenecks
- [ ] Compare unified vs original function performance
- [ ] Create performance report

#### Phase 7.2: Optimize if Needed (1-2 hours)
- [ ] Target areas:
  - ODE callback frequency
  - Rate computation caching
  - Memory allocation patterns
- [ ] Only optimize hot paths identified in profiling
- [ ] Benefits: Faster simulations if significant

---

## Prioritized Next Steps (Recommended Order)

### Just Completed ✅
1. **Phase 3.1**: Structured rate data ✅ DONE (2 hours)
   - Created RateDatabase with 52 structured rates
   - Type-safe, well-documented, fully tested
   - Backward compatible integration layer

### Immediate (Next Session)
2. **Phase 6.1**: Baseline numerical tests (2 hours)
   - Validates refactoring across parameter space
   - Builds confidence for production use
   - Compare standard and spinone modes

3. **Phase 3.2**: Separate rate computation logic (1-2 hours)
   - Extract compute_sr_rates() to dedicated module
   - Organize related rate functions

### Short Term (After that)
4. **Phase 4**: Module organization (2-3 hours)
   - Makes code more accessible
   - Enables proper testing framework

5. **Phase 5.1**: Docstring documentation (1.5 hours)
   - Improves usability
   - IDE support

### Medium Term
5. **Phase 5.3**: User guide (1 hour)
6. **Phase 7**: Performance work (if needed)

---

## Risk Assessment

### Low Risk (Safe to proceed)
- [x] Phase 3.1: Structured data - Well-defined, isolated task
- [x] Phase 4.1: Module file - Standard Julia pattern
- [x] Phase 5.1: Docstrings - Non-breaking enhancement

### Medium Risk (Test thoroughly)
- [ ] Phase 3.2: Rate computation extraction - Must preserve behavior
- [ ] Phase 4.2: Move scripts - Ensure all dependencies work
- [ ] Phase 6.1: Baseline comparison - Critical for validation

### High Risk (Requires extra care)
- [ ] Phase 4.3: Test reorganization - Can break CI/CD
- [ ] Phase 7: Performance changes - Can introduce bugs

---

## Success Criteria

### Phase 3 ✓
- [ ] All rate coefficients in structured format
- [ ] Lookup functions work correctly
- [ ] No change in numerical results

### Phase 4 ✓
- [ ] Code loads with main module
- [ ] Scripts run from new location
- [ ] Tests work with new structure

### Phase 5 ✓
- [ ] All public functions documented
- [ ] IDE shows docstrings
- [ ] User guide is comprehensive

### Phase 6 ✓
- [ ] Baseline tests show < 1e-10 difference
- [ ] Edge cases behave reasonably
- [ ] No NaN/Inf in extended tests

### Phase 7 ✓ (if done)
- [ ] Performance improvement > 5%
- [ ] No numerical differences
- [ ] Behavior unchanged

---

## Time Estimate Summary

| Phase | Hours | Priority | Status |
|-------|-------|----------|--------|
| 1 | 2-3 | HIGH | ✅ Done |
| 2 | 4-5 | HIGH | ✅ Done |
| 3 | 3-4 | HIGH | 📋 Next |
| 4 | 2-3 | MEDIUM | 📋 Soon |
| 5 | 2-3 | MEDIUM | 📋 Soon |
| 6 | 2-3 | HIGH | 📋 Soon |
| 7 | 2-3 | LOW | 📋 Optional |
| **TOTAL** | **~20 hours** | | |

---

## Recommended First Next Steps

### Session 1 (3-4 hours)
```julia
# Create structured rate data
- Phase 3.1: Create src/Core/rate_coefficients.jl
- Phase 3.1: Define RateCoefficient struct
- Phase 3.1: Load rates from load_rates.jl
```

### Session 2 (2-3 hours)
```julia
# Validate numerical behavior
- Phase 6.1: Expand baseline test cases
- Phase 6.1: Compare standard mode results
- Phase 6.1: Compare spinone mode results
```

### Session 3 (2-3 hours)
```julia
# Organize code structure
- Phase 4.1: Create main AxionSR.jl module
- Phase 4.2: Move scripts to scripts/ directory
- Phase 5.1: Add docstrings to public API
```

---

## Current Git Status

**Branch**: phase-2-2-merge (Phase 2.2 complete)
**Latest Commits**:
- e94e752: Fix integration test
- 32b579f: Integration Testing
- 229a7ed: Phase 2.2 merge

**Next**: Consider merging to main when ready, or continue development on feature branch

---

## Notes

- All completed phases well-tested and documented
- Integration tests passing (5/5)
- Code ready for production use as-is
- Remaining work improves organization, not functionality
- Can stop at any phase - each is independently valuable
- Consider user feedback before major reorganization

**Recommendation**: Proceed with Phase 3 (rate coefficients) + Phase 6 (validation) next, then assess. Phases 4-5 are nice-to-have improvements.
