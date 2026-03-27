# AxionSR User Manual - Complete Index

## 📚 All Files at a Glance

### Core Documentation Files

| File | Type | Lines | Purpose |
|------|------|-------|---------|
| **main.tex** | LaTeX | 320 | Master document with preamble |
| **references.bib** | BibTeX | 200 | Complete bibliography (50+ refs) |

### Main Content Sections

| File | Lines | Topics |
|------|-------|--------|
| **sections/01_theory.tex** | 280 | Axion physics, superradiance, quantum evolution |
| **sections/02_methods.tex** | 410 | Eigenvalue solvers, rate computation, ODE integration, MCMC |
| **sections/03_usage.tex** | 360 | Installation, API reference, workflows, examples |
| **sections/04_examples.tex** | 350 | 5 detailed worked examples with code |

### Appendices

| File | Lines | Topics |
|------|-------|--------|
| **appendices/A_mathematics.tex** | 200 | Mathematical derivations, proofs, equations |
| **appendices/B_rate_tables.tex** | 200 | Data format, file structure, interpolation |
| **appendices/C_troubleshooting.tex** | 250 | FAQ, debugging, optimization, self-checks |

### Quick Reference & Guides

| File | Type | Purpose |
|------|------|---------|
| **README.md** | Markdown | Manual overview, build instructions |
| **QUICK_REFERENCE.md** | Markdown | One-page cheat sheet (printable) |
| **MANUAL_SUMMARY.md** | Markdown | Complete summary and statistics |
| **INDEX.md** | Markdown | This file - navigation guide |

## 🎯 Finding What You Need

### By User Type

#### **First-time Users**
1. Read: **README.md** (2 min)
2. Read: **Section 1: Theory** (15 min)
3. Read: **Section 3.1: Installation** (5 min)
4. Run: **Example 1** (10 min)
5. Print: **QUICK_REFERENCE.md**

#### **Implementing Workflows**
1. Reference: **Section 3.3: Common Workflows** (10 min)
2. Study: **Example 2-3** (30 min)
3. Adapt code for your system (60 min)
4. Consult: **Appendix C** if needed

#### **Researchers Deriving Constraints**
1. Study: **Section 2: Methods** (30 min)
2. Study: **Example 3: MCMC Inference** (30 min)
3. Review: **Section 1.4: Observables** (15 min)
4. Reference: **Appendix A: Mathematics** as needed
5. Cite: **references.bib** for publications

#### **Developers & Contributors**
1. Study: **Section 2: Methods** (60 min)
2. Review: **Appendix B: Rate Tables** (30 min)
3. Read: **Source code** in main repository
4. Reference: **Appendix A: Mathematics** for validation
5. Check: **Section 3.4: Advanced API**

#### **Debugging & Troubleshooting**
1. First check: **QUICK_REFERENCE.md** - Troubleshooting table
2. Reference: **Appendix C.1-C.4** - Detailed diagnostics
3. Search: **Section 3: Usage** for API details
4. Check: **Appendix C.5: FAQs**
5. Report: Issue to GitHub with diagnostics from C.4

### By Topic

#### Axion Physics
- Superradiance condition: **Section 1.1**
- Teukolsky equation: **Section 1.2**
- Quantum evolution: **Section 1.3**
- Observable signatures: **Section 1.4**
- Mathematical details: **Appendix A.1-A.6**

#### Computational Methods
- Eigenvalue solver: **Section 2.1** + **Appendix A.3-A.4**
- Rate computation: **Section 2.2** + **Appendix B**
- Time integration: **Section 2.3**
- Statistical inference: **Section 2.4**
- Performance: **Section 2.5**

#### Implementation
- Installation: **Section 3.1**
- Core functions: **Section 3.2**
- Workflows: **Section 3.3** + **Section 4**
- Advanced features: **Section 3.4**
- Configuration: **Section 3.5**

#### Examples & Data
- Cygnus X-1: **Example 1** (Section 4.1)
- Population analysis: **Example 2** (Section 4.2)
- MCMC inference: **Example 3** (Section 4.3)
- Bosenova signatures: **Example 4** (Section 4.4)
- Sensitivity analysis: **Example 5** (Section 4.5)

#### Troubleshooting
- Installation: **Appendix C.1**
- Computation errors: **Appendix C.2**
- Performance: **Appendix C.3**
- Verification: **Appendix C.4**
- Common questions: **Appendix C.5**
- Help resources: **Appendix C.6**

## 📖 Reading Paths by Goal

### **Goal: Learn Axion Superradiance**
```
README → Section 1 → QUICK_REFERENCE → Appendix A
Time: 2 hours
```

### **Goal: Use Code for My Black Hole System**
```
Section 3.1 → Section 3.2 → Example 1 → Example 2 → Modify for your data
Time: 3 hours
```

### **Goal: Derive Constraints on Axion Mass**
```
Section 1 → Section 2.4 → Example 2 → Example 3 → Advanced MCMC
Time: 6 hours
```

### **Goal: Understand the Physics in Depth**
```
Section 1 → Appendix A → Section 2 → Appendix B → Section 4
Time: 10+ hours
```

### **Goal: Fix an Error or Unexpected Result**
```
QUICK_REFERENCE (Troubleshooting) → Appendix C.1-C.4 → Section 3
Time: 30 minutes to 2 hours
```

## 🔗 Cross-Reference Map

```
main.tex
├─→ Section 1: Theory
│   ├─→ Section 1.1: Axion superradiance → Appendix A.1-A.2
│   ├─→ Section 1.2: Teukolsky equation → Appendix A.3-A.4
│   ├─→ Section 1.3: Quantum evolution → Appendix A.5-A.6
│   └─→ Section 1.4: Observables → Example 2-3
│
├─→ Section 2: Methods
│   ├─→ Section 2.1: Eigenvalue → Appendix A.3-A.4
│   ├─→ Section 2.2: Rate computation → Appendix B
│   ├─→ Section 2.3: Time integration → Appendix A.5
│   ├─→ Section 2.4: Statistical inference → Appendix A.8
│   └─→ Section 2.5: Performance → Section 3.4
│
├─→ Section 3: Usage Guide
│   ├─→ Section 3.1: Installation → README.md
│   ├─→ Section 3.2: API → QUICK_REFERENCE.md
│   ├─→ Section 3.3: Workflows → Section 4
│   ├─→ Section 3.4: Advanced → Section 2
│   └─→ Section 3.5: Configuration → Appendix B
│
├─→ Section 4: Examples
│   ├─→ Example 1: Single BH → Section 1.4
│   ├─→ Example 2: Population → Section 1.4, Section 2.4
│   ├─→ Example 3: MCMC → Section 2.4, Appendix A.8
│   ├─→ Example 4: Bosenova → Section 1.3
│   └─→ Example 5: Sensitivity → Section 1, Section 2
│
└─→ Appendices
    ├─→ Appendix A: Mathematics → Referenced throughout
    ├─→ Appendix B: Rate Tables → Section 2.2, Section 3.5
    └─→ Appendix C: Troubleshooting → All sections
```

## 📝 How to Use This Manual

### **Reading the PDF (After Compilation)**
- Use PDF viewer bookmark feature for quick navigation
- Click hyperlinks to jump between sections
- Use PDF search to find topics
- Print QUICK_REFERENCE.md for quick lookup while coding

### **Citing Specific Sections**
- Equations: `Eq. (1.23)` or `Equation 1.23`
- Sections: `Section 2.1` or `§2.1`
- Examples: `Example 3` or `Section 4.3`
- Appendices: `Appendix A.4` or `§A.4`

### **Contributing to Manual**
To submit improvements:
1. Edit relevant `.tex` file
2. Rebuild: `latexmk -pdf main.tex`
3. Submit pull request with rationale
4. Update version in main.tex preamble

## 📊 Manual Statistics

| Metric | Value |
|--------|-------|
| **Total files** | 12 |
| **Total lines** | 2,751 |
| **Total size** | 136 KB |
| **Main sections** | 4 |
| **Appendices** | 3 |
| **Code examples** | 20+ |
| **Figures/tables** | 15+ |
| **Bibliography entries** | 50+ |
| **Estimated PDF pages** | 60-80 |

## 🔐 File Locations

All files located in: `/Users/samuelwitte/Dropbox/Axion_SR/UserManual/`

```
UserManual/
├── main.tex                    ← Start here
├── references.bib
├── README.md
├── QUICK_REFERENCE.md
├── MANUAL_SUMMARY.md
├── INDEX.md (this file)
├── sections/
│   ├── 01_theory.tex
│   ├── 02_methods.tex
│   ├── 03_usage.tex
│   └── 04_examples.tex
└── appendices/
    ├── A_mathematics.tex
    ├── B_rate_tables.tex
    └── C_troubleshooting.tex
```

## ✅ Verification Checklist

Before using manual, verify:
- [ ] All `.tex` files present (10 files)
- [ ] All `.bib` file present (1 file)
- [ ] All `.md` guides present (4 files)
- [ ] Directory structure matches above
- [ ] Can compile: `latexmk -pdf main.tex`
- [ ] PDF generates without errors
- [ ] Table of contents accurate
- [ ] All cross-references valid

## 🎓 Suggested Use in Papers

### In Methodology Sections
> "We use the AxionSR computational framework (Witte 2025) to simulate axion superradiance evolution..."

### As Supplementary Material
> "Technical details of eigenvalue computation and rate interpolation are provided in the AxionSR User Manual (Witte 2025), available at [URL]"

### For Data & Methods
Reference specific sections:
- Section 1: For physics background
- Section 2: For computational method details
- Section 3: For reproducibility
- Section 4: For validation examples

## 📞 Support & Resources

| Resource | Location | Use |
|----------|----------|-----|
| **Quick help** | QUICK_REFERENCE.md | Fast answers |
| **Installation** | Section 3.1 | Getting started |
| **API docs** | Section 3.2 | Function reference |
| **Examples** | Section 4 | Learn by doing |
| **Theory** | Section 1 | Background physics |
| **Troubleshooting** | Appendix C | Debugging |
| **Advanced topics** | Appendix A-B | Deep dive |
| **Bibliography** | references.bib | Literature |

---

**Last Updated**: November 23, 2025
**Version**: 1.0
**Status**: ✅ Complete and ready for use
