# AxionSR Open-Source Paper - Quick Start Guide

## Files Created

- **`paper.tex`** - Full arXiv-ready LaTeX paper (~500 lines)
- **`OVERLEAF_SETUP.md`** - Detailed Overleaf setup instructions
- **This file** - Quick reference

## What's in the Paper

The template includes 10+ sections covering:

| Section | Status | Notes |
|---------|--------|-------|
| Title & Abstract | ✅ Done | Placeholder content - customize with your details |
| Introduction | ✅ Structure | Key motivation, related work references |
| Physical Model | ✅ Structure | Superradiance, dynamics, bosenova |
| Numerical Methods | ✅ Structure | Eigenvalue, rates, time evolution, boundary conditions |
| Code Structure | ✅ Structure | File organization, main functions, architecture |
| Validation | ✅ Structure | Analytical limits, benchmarks, tests |
| Features | ✅ Structure | Physical modes, configuration options |
| Performance | ✅ Structure | Runtime, accuracy, stability |
| Applications | ✅ Structure | Dark matter, black hole demographics, extensions |
| Software Engineering | ✅ Structure | Code quality, reproducibility, open science |
| Installation & Usage | ✅ Structure | Quick start, examples |
| Future Directions | ✅ Structure | Planned enhancements |
| Appendices | ✅ Structure | Math formalism, examples |

## Quickest Path to Overleaf

### Step 1: Create Overleaf Account
Go to https://www.overleaf.com/ and sign up (if needed)

### Step 2: Create New Project
- Click **New Project** → **Blank Project**
- Name: "AxionSR - Axion Superradiance Code"

### Step 3: Upload Paper
**Option A - Copy-Paste (2 minutes)**:
1. Open `paper.tex` locally
2. Select all (Ctrl+A)
3. Copy (Ctrl+C)
4. In Overleaf, delete default text
5. Paste (Ctrl+V)
6. Hit green **Recompile** button

**Option B - GitHub Sync (5 minutes)**:
1. Push your code to GitHub
2. In Overleaf: **New Project** → **Import from GitHub**
3. Select `YOUR_USERNAME/AxionSR`
4. Auto-syncs from now on!

### Step 4: Customize
Find and replace these in `paper.tex`:

```
FIND                        REPLACE WITH
========================================
Samuel Witte                Your Name
samuelwitte@physics.ucla.edu Your Email
[collaborators/advisors]     Actual names
[funding sources]            Your funding
```

### Step 5: Add Content
- Expand Methods section with your implementation details
- Add Validation & Results sections with plots
- Update References with actual citations
- Insert figures (in Overleaf: **Files** → **Upload**)

---

## Key Sections to Personalize

### 1. Title Page (Lines 40-48)
```latex
\title{AxionSR: An Open-Source Code for Axion Superradiance Dynamics}
\author{Samuel Witte\thanks{...}}
```
→ Replace with your names and affiliations

### 2. Abstract (Lines 51-62)
Currently: Generic template
→ Write your specific abstract (250-300 words)

### 3. Introduction (Lines 65-120)
- Motivation ✅ (Structure present)
- Challenge ✅ (Structure present)
- Your contribution → ADD YOUR DETAILS

### 4. Methods (Lines 123-220)
- Physical model ✅ (Equations present)
- Numerical methods ✅ (Algorithms described)
- Implementation details → EXPAND WITH YOUR CODE

### 5. New Sections to Add

**Results Section** (after Validation):
```latex
\section{Results and Benchmarks}

\subsection{Superradiance Timescales}
% Add figure and discussion

\subsection{Population Studies}
% Add your parameter space exploration results
```

**Figures Section** (in Methods):
```latex
\begin{figure}[h]
    \centering
    \includegraphics[width=0.8\textwidth]{path/to/your/figure.pdf}
    \caption{Your figure caption.}
    \label{fig:your_label}
\end{figure}
```

### 6. References (Lines 451-475)
Currently: Placeholder citations
→ Replace with actual references using this format:

```latex
\bibitem{FirstAuthor2024} First Author, et al. (2024).
Paper Title. \textit{Journal Name}, Vol(Issue), Pages.
arXiv:XXXX.XXXXX.
```

### 7. Acknowledgments (Lines 448-449)
```latex
\section*{Acknowledgments}

We thank [names] for helpful discussions. This work was
supported by [funding agencies/grants].
```

---

## Adding Figures to Overleaf

### Method 1: Upload Files
1. Click **Files** (left sidebar)
2. Click **Upload** button
3. Select your PDF/PNG files
4. Place in document with:
```latex
\begin{figure}[h]
    \centering
    \includegraphics[width=0.8\textwidth]{filename.pdf}
    \caption{Caption here.}
    \label{fig:label}
\end{figure}
```

### Method 2: Create Inline
For simple diagrams, use TikZ (built into Overleaf):
```latex
\begin{tikzpicture}
  \draw (0,0) -- (1,0) -- (1,1) -- (0,1) -- cycle;
\end{tikzpicture}
```

---

## Common LaTeX Patterns in Your Paper

### Equations
```latex
% Single equation (numbered)
\begin{equation}
\alpha = \mu M G_N
\end{equation}

% Inline equation
The parameter $\alpha$ controls superradiance efficiency.

% Multiple equations
\begin{align}
E &= mc^2 \\
p &= \gamma m v
\end{align}
```

### Citations
```latex
% Numerical citation
As shown in \cite{Arvanitaki2015}...

% Multiple citations
Several works \cite{Brito2015,East2017} have explored...
```

### Cross-References
```latex
% Reference to section
See Section~\ref{sec:intro} for background.

% Reference to equation
From Equation~\eqref{eq:mass_evolution}...

% Reference to figure
Figure~\ref{fig:spin_evolution} shows...
```

### Code Listings
Already configured in preamble! Use:
```latex
\begin{lstlisting}[language=Julia, caption=Your Caption]
function solve_system(mu, fa, aBH, M_BH, t_max)
    # Your code here
end
\end{lstlisting}
```

### Tables
```latex
\begin{table}[h]
    \centering
    \caption{Performance Metrics}
    \begin{tabular}{lcc}
        \toprule
        Method & Time (s) & Accuracy \\
        \midrule
        Heun & 0.5 & $10^{-8}$ \\
        Cheby & 0.3 & $10^{-10}$ \\
        \bottomrule
    \end{tabular}
    \label{tab:performance}
\end{table}
```

---

## Before Submitting to arXiv

### Checklist

- [ ] Title and authors finalized
- [ ] Abstract is clear and complete (250-300 words)
- [ ] Introduction motivates the problem
- [ ] Methods section is detailed and clear
- [ ] Results section with figures and analysis
- [ ] Validation against known results
- [ ] All citations are complete and accurate
- [ ] Equations are numbered correctly
- [ ] Figures have captions and labels
- [ ] References section complete
- [ ] PDF compiles without errors
- [ ] Grammar and spelling checked
- [ ] Formatting is consistent (fonts, spacing)

### Export from Overleaf

1. Click **Menu** (top left)
2. Download → **PDF**
3. Or: Source files (for arXiv submission)
4. Or: Git clone (for version control)

---

## arXiv Submission Process

### 1. Prepare Files
```bash
# Create a folder with all files
mkdir arxiv_submission
cd arxiv_submission

# Copy main paper
cp ../paper.tex .

# Copy all figures (must be in same directory for arXiv)
cp ../plts/*.pdf .

# Create .bbl file from bibliography (Overleaf downloads this)
# Download from Overleaf menu → Source
```

### 2. Go to arXiv
- Visit https://arxiv.org/
- New authors: **Register**
- Submit article

### 3. Fill in Metadata
- **Title**: "AxionSR: An Open-Source Code for..."
- **Authors**: Your name(s)
- **Abstract**: Copy from paper
- **Categories**:
  - Primary: gr-qc (General Relativity)
  - Secondary: astro-ph.HE (High-Energy Astrophysics)
- **Subjects**: Black holes, Axions, Numerical Relativity

### 4. Upload Files
- Paper.tex
- All .pdf figures
- .bbl file (if needed)

### 5. Review and Submit
- Check preview PDF
- Click Submit
- Get your arXiv ID!

---

## Helpful Overleaf Keyboard Shortcuts

| Action | Shortcut |
|--------|----------|
| Recompile | Ctrl+Enter |
| Find | Ctrl+F |
| Find & Replace | Ctrl+H |
| Bold | Ctrl+B |
| Italic | Ctrl+I |
| Comment out line | Ctrl+/ |

---

## Getting Help

### Overleaf Support
- **Overleaf Help**: https://www.overleaf.com/help
- **Learn LaTeX**: https://www.overleaf.com/learn
- **Chat with support**: Help button (?) in Overleaf

### LaTeX Resources
- **Stack Exchange**: https://tex.stackexchange.com/
- **r/LaTeX**: https://reddit.com/r/LaTeX
- **CTAN**: https://www.ctan.org/ (packages)

### Physics/Math Resources
- **arXiv Help**: https://arxiv.org/help/
- **INSPIRE HEP**: https://inspirehep.net/ (citations)
- **NASA ADS**: https://ui.adsabs.harvard.edu/ (papers)

---

## Next Steps

**Right now**:
1. Go to https://www.overleaf.com/
2. Create account or log in
3. Copy `paper.tex` content into a new project
4. Click Recompile and see it work!

**This week**:
1. Customize title, authors, abstract
2. Expand Methods section
3. Add your Results

**When ready**:
1. Polish and peer review
2. Submit to arXiv
3. Get your arXiv ID!
4. Share with community!

---

Enjoy writing your paper! 🚀
