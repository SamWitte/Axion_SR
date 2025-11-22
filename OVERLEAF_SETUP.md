# Setting Up Your Paper on Overleaf

This guide will help you upload and sync your AxionSR paper with Overleaf for collaborative editing.

## Option 1: Upload LaTeX File Directly (Quickest)

1. **Go to Overleaf**: https://www.overleaf.com/
2. **Create Account** (if needed) and log in
3. **New Project** → **Blank Project**
4. **Project Name**: `AxionSR - Axion Superradiance Code`
5. **Copy & Paste Method**:
   - Open `paper.tex` in your text editor
   - Copy all content (Ctrl+A, Ctrl+C)
   - In Overleaf, delete the default content
   - Paste the full paper.tex (Ctrl+V)
   - Click **Recompile** (green button)

✅ **Pros**: Quick, works immediately
❌ **Cons**: Changes not synced to GitHub

---

## Option 2: GitHub-Connected Project (Recommended)

This allows two-way sync between your GitHub repo and Overleaf.

### Step 1: Set Up GitHub Repository

If you don't have one already:

```bash
cd /Users/samuelwitte/Dropbox/Axion_SR
git remote add origin https://github.com/YOUR_USERNAME/AxionSR.git
git branch -M main
git push -u origin main
```

### Step 2: Connect Overleaf to GitHub

1. Go to **https://www.overleaf.com/home**
2. Click **New Project** → **Import from GitHub**
3. You'll be prompted to authorize Overleaf to access your GitHub
4. Click **Authorize Overleaf**
5. Select `YOUR_USERNAME/AxionSR` from the dropdown
6. Click **Import**

### Step 3: Two-Way Syncing

**From GitHub to Overleaf** (automatic):
- Push changes to GitHub: `git push`
- Overleaf will sync automatically

**From Overleaf to GitHub**:
1. Make edits in Overleaf
2. Click **Menu** (top left) → **GitHub** → **Push overleaf/main to GitHub**
3. Create a pull request or commit directly

---

## Option 3: Git + SSH (Advanced)

For automated syncing without the GitHub middleman:

1. In Overleaf, click **Menu** → **Sync** → **Git**
2. Follow instructions to set up SSH keys
3. Clone to local machine:
```bash
git clone git@git.overleaf.com:PROJECT_ID/PROJECT_NAME.git overleaf-sync
```

---

## Editing Your Paper

### Key Sections to Update

Open `paper.tex` and find these sections to customize:

1. **Title & Authors** (line ~40)
   ```latex
   \title{AxionSR: An Open-Source Code for Axion Superradiance Dynamics}
   \author{Your Name, Your Affiliation}
   ```

2. **Abstract** (line ~55)
   - Replace placeholder with your own abstract

3. **Introduction** (line ~85)
   - Add your motivation and related work
   - Update citations as needed

4. **Methods** (line ~168)
   - Expand sections with implementation details
   - Add equations, algorithms, figures

5. **Results** (new section)
   - Add validation results
   - Include figures comparing to benchmarks
   - Performance metrics

6. **Acknowledgments** (line ~445)
   - Thank collaborators, advisors, funding agencies

7. **References** (line ~447)
   - Update bibliography with your actual citations
   - Use BibTeX format or natbib

### Adding Figures

To insert a figure after results section:

```latex
\begin{figure}[h]
    \centering
    \includegraphics[width=0.8\textwidth]{figs/superradiance_evolution.pdf}
    \caption{Black hole spin evolution during superradiance.}
    \label{fig:spin_evolution}
\end{figure}
```

Then reference it with: `As shown in Fig.~\ref{fig:spin_evolution}...`

### Adding Citations

Update the bibliography section with proper references:

```latex
\bibitem{Your2024} Your Name, et al. (2024). Your Paper Title.
\textit{Journal Name}, Volume(Issue), pages. arXiv:xxxx.xxxxx.
```

Then cite with: `\cite{Your2024}`

---

## Compiling & Preview

1. **Recompile**: Click green "Recompile" button (top right)
2. **View PDF**: Preview pane on right side
3. **Full Screen PDF**: Click expand icon

### Troubleshooting Compilation

- **Missing packages**: Overleaf has most common packages
- **Bibliography errors**: Ensure BibTeX format is correct
- **Figure not found**: Check file names match exactly
- **Encoding issues**: Add `\usepackage[utf-8]{inputenc}` (already included)

---

## Customization Options

### Fonts
Current: Computer Modern (default LaTeX)

To use different fonts, add after `\documentclass`:
```latex
\usepackage{times}      % Times New Roman
\usepackage{mathptmx}   % Math in Times
```

### Page Layout
Current: 1-inch margins, 11pt font

Change with:
```latex
\geometry{margin=0.75in}  % narrower margins
\documentclass[12pt]{article}  % larger font
```

### Color Scheme
Current: Blue hyperlinks

Change link colors:
```latex
\hypersetup{
    colorlinks=true,
    linkcolor=darkblue,
    citecolor=darkgreen,
    urlcolor=darkred
}
```

---

## Submitting to arXiv

Once your paper is ready:

1. **Export as PDF** from Overleaf: Menu → Download → PDF
2. **Create arXiv account**: https://arxiv.org/
3. **Submit paper**:
   - Title, abstract, authors
   - Upload source files (LaTeX + figures)
   - OR upload compiled PDF
4. **Choose category**:
   - astro-ph.HE (high-energy astrophysics)
   - gr-qc (general relativity)
   - hep-ph (high-energy physics - phenomenology)
5. **Submit and receive arXiv ID**

Your paper will be assigned an ID like: `2410.12345`

---

## Collaboration Features

### Invite Collaborators

In Overleaf:
1. Click **Share** (top right)
2. **Edit link**: Can anyone with link edit it
3. **Email**: Add specific collaborators by email
4. **Permissions**: Choose view-only or edit

### Track Changes

Overleaf has built-in:
- **Edit history**: View past versions
- **Comments**: Add reviewer comments
- **Change tracking**: See who changed what when

### Keeping GitHub in Sync

Best workflow:
```bash
# Make local edits
# Overleaf syncs to GitHub automatically (if connected)
# Or manually push from Overleaf

# Pull latest changes
git pull origin main

# Continue working
```

---

## Tips for Writing Physics Papers

1. **Equation formatting**: Use `\begin{equation}` for important equations, `\(` `\)` for inline math
2. **References**: Use `\label{sec:something}` and `\ref{sec:something}` for cross-references
3. **Subsections**: Organize with `\section`, `\subsection`, `\subsubsection`
4. **Tables**: Use `\begin{table}` environment with `\begin{tabular}`
5. **Code listings**: Use `\begin{lstlisting}` for code samples (already configured)

---

## Useful Overleaf Resources

- **Overleaf Learn**: https://www.overleaf.com/learn
- **LaTeX Templates**: https://www.overleaf.com/templates
- **Documentation**: https://www.overleaf.com/learn/latex/
- **Quick Reference**: https://www.overleaf.com/learn/latex/Learn_LaTeX_in_30_minutes

---

## Next Steps

1. ✅ Copy `paper.tex` to Overleaf (or connect via GitHub)
2. ✅ Customize title, authors, abstract
3. ✅ Expand methods and results sections
4. ✅ Add your figures and results
5. ✅ Update citations and references
6. ✅ Compile and review PDF
7. ✅ Invite collaborators for feedback
8. ✅ Submit to arXiv!

---

**Questions?** Check Overleaf's help menu or visit the Overleaf community forums.
