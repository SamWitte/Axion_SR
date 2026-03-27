# AxionSR User Manual

Professional user manual and technical guide for the **AxionSR** computational framework for simulating axion superradiance in rotating black holes.

## Contents

This manual provides a comprehensive reference for researchers and practitioners using AxionSR to:
- Compute eigenfrequencies of quasi-bound states in Kerr spacetime
- Simulate long-term evolution of axion condensates around spinning black holes
- Derive astrophysical constraints on ultralight boson properties
- Analyze multi-messenger observations from black hole populations

## Document Structure

### Main Sections

1. **Introduction** (`main.tex`)
   - Physical motivation for axion superradiance
   - Scope and organization of the code
   - Manual overview

2. **Section 1: Theoretical Framework** (`sections/01_theory.tex`)
   - Axion physics and superradiance condition
   - Teukolsky equation and quasi-bound states
   - Quantum state evolution and nonlinear effects
   - Astrophysical observables and constraints

3. **Section 2: Computational Methods** (`sections/02_methods.tex`)
   - Eigenvalue computation via Heun equation
   - Rate coefficient computation
   - Time integration with callbacks
   - Statistical inference (MCMC)
   - Performance analysis

4. **Section 3: User Guide and API** (`sections/03_usage.tex`)
   - Installation instructions
   - Core API reference
   - Common workflows
   - Advanced usage patterns
   - Performance tips and optimization

5. **Section 4: Worked Examples** (`sections/04_examples.tex`)
   - Example 1: Cygnus X-1 spin prediction
   - Example 2: X-ray binary population constraints
   - Example 3: Bayesian parameter inference
   - Example 4: Bosenova gravitational wave signatures
   - Example 5: Sensitivity analysis

### Appendices

- **Appendix A: Mathematical Details** (`appendices/A_mathematics.tex`)
  - Teukolsky equation derivation
  - Heun recursion relation
  - Boundary conditions and eigenvalue computation
  - Bayesian posterior formulas

- **Appendix B: Rate Table Structure** (`appendices/B_rate_tables.tex`)
  - Pre-computed rate file format
  - NPZ file structure
  - Interpolation procedures
  - Custom rate table generation

- **Appendix C: Troubleshooting and FAQs** (`appendices/C_troubleshooting.tex`)
  - Installation issues
  - Computation failures
  - Performance optimization
  - Results verification
  - Common questions

### References

- **Bibliography** (`references.bib`)
  - Axion physics literature
  - Superradiance theory
  - Black hole astrophysics
  - Numerical methods
  - Gravitational wave astronomy

## Building the Manual

### Prerequisites

- **LaTeX distribution**: TeX Live (Linux/macOS) or MiKTeX (Windows)
- **Required packages**:
  - geometry, amsmath, amssymb, listings, hyperref, natbib, float, algorithm, algorithmicx

### Compilation

```bash
# One-shot compilation
pdflatex main.tex
bibtex main.aux
pdflatex main.tex
pdflatex main.tex

# Or using latexmk (recommended)
latexmk -pdf main.tex
```

### Output

The compiled PDF will be generated as `main.pdf` (typically 50-100 pages).

## File Organization

```
UserManual/
├── main.tex                     # Main document file
├── sections/
│   ├── 01_theory.tex           # Theoretical foundations
│   ├── 02_methods.tex          # Computational methods
│   ├── 03_usage.tex            # User guide and API
│   └── 04_examples.tex         # Worked examples
├── appendices/
│   ├── A_mathematics.tex       # Mathematical details
│   ├── B_rate_tables.tex       # Rate table documentation
│   └── C_troubleshooting.tex   # Troubleshooting guide
├── references.bib              # Bibliography
├── README.md                   # This file
└── main.pdf                    # Compiled output (after building)
```

## Usage

This manual is intended for:

1. **New users**: Start with Section 1 and Section 3 (Usage Guide)
2. **Developers**: Refer to Section 2 (Methods) and Appendices for implementation details
3. **Researchers**: Use Sections 4 and references for application examples
4. **Troubleshooting**: Consult Appendix C and the main README in the code repository

## Citation

Please cite this manual in academic work using:

```bibtex
@manual{AxionSR2025,
  author = {Witte, Samuel D.},
  title = {AxionSR: A Computational Framework for Axion Superradiance
           in Black Holes -- User Manual and Technical Guide},
  year = {2025},
  url = {https://github.com/your-repo/AxionSR.jl}
}
```

## License

This manual is provided under the same license as the AxionSR codebase.

## Contributing

Suggestions for improvements to this manual are welcome. Please submit issues or pull requests to the main repository.

## Contact

For questions about the manual or code:
- **GitHub Issues**: [AxionSR Issues](https://github.com/your-repo/AxionSR.jl/issues)
- **Email**: [Contact information]

---

**Last Updated**: November 2025
**Manual Version**: 1.0
**Corresponding Code Version**: AxionSR v1.0+
