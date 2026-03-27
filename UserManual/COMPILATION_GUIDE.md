# Compilation Guide for AxionSR User Manual

## Prerequisites

### macOS
```bash
# Install MacTeX (includes pdflatex, bibtex, latexmk)
brew install mactex
# Or download from: https://tug.org/mactex/

# Verify installation
pdflatex --version
latexmk -version
```

### Linux (Ubuntu/Debian)
```bash
# Install TeX Live
sudo apt-get install texlive-latex-full texlive-fonts-recommended texlive-xetex latexmk

# Verify
pdflatex --version
latexmk -version
```

### Windows
- Download MiKTeX: https://miktex.org/download
- Or download TeX Live: https://tug.org/texlive/
- Both include all necessary tools

## Quick Compilation

### Method 1: Using latexmk (Recommended)

The simplest and most reliable method:

```bash
cd /Users/samuelwitte/Dropbox/Axion_SR/UserManual
latexmk -pdf main.tex
```

This automatically:
- Runs pdflatex multiple times
- Processes bibliography (bibtex)
- Updates cross-references
- Generates final PDF

**Expected output:**
```
Latexmk: applying rule 'pdflatex'...
...
pdflatex: successfully generated 'main.pdf'
Latexmk: Finished successfully
```

The PDF should appear as `main.pdf` in the same directory.

### Method 2: Manual Step-by-Step

If latexmk is unavailable:

```bash
# Step 1: First LaTeX pass
pdflatex main.tex

# Step 2: Process bibliography
bibtex main.aux

# Step 3: Second LaTeX pass (update references)
pdflatex main.tex

# Step 4: Third LaTeX pass (finalize cross-references)
pdflatex main.tex
```

After step 4, `main.pdf` will be generated.

## Troubleshooting Compilation

### Problem: "pdflatex: command not found"

**Solution**: LaTeX is not installed
- macOS: Run `brew install mactex`
- Linux: Run `sudo apt-get install texlive-latex-full`
- Windows: Download and install MiKTeX or TeX Live

### Problem: "Package not found: siunitx" (or other packages)

**Solution**: Missing LaTeX packages
```bash
# macOS (with MacTeX)
# Run TeX Live Utility and install missing packages
# Or use: tlmgr install siunitx

# Linux
sudo apt-get install texlive-fonts-recommended

# Windows (MiKTeX)
# Run MiKTeX Console and install missing packages
```

### Problem: Bibliography doesn't appear

**Ensure you run all 4 steps** (or use latexmk which does this automatically)
- Step 1: `pdflatex main.tex`
- Step 2: `bibtex main.aux`
- Step 3: `pdflatex main.tex`
- Step 4: `pdflatex main.tex`

### Problem: "File not found: sections/01_theory.tex"

**Solution**: Run from correct directory
```bash
# Wrong:
cd ~ && pdflatex /path/to/UserManual/main.tex

# Correct:
cd /Users/samuelwitte/Dropbox/Axion_SR/UserManual
pdflatex main.tex
```

### Problem: Cross-references show "??" or "p. ?"

**Solution**: Run LaTeX multiple times
```bash
# Run a third pass:
pdflatex main.tex
```

## Verification

After compilation, verify the PDF is correct:

```bash
# Check file was created and is reasonable size
ls -lh main.pdf
# Should show ~500 KB to 1 MB

# Check it's a valid PDF
file main.pdf
# Should show: "PDF document, version 1.5"

# Verify page count
pdfinfo main.pdf | grep Pages
# Should show ~60-80 pages
```

## Advanced Options

### Generate bookmarks/hyperlinks

The manual includes hyperlinks by default (hyperref package).
Open main.pdf in any PDF viewer and use the bookmarks panel.

### Change paper size or margins

Edit `main.tex` line 5:
```latex
\usepackage[margin=1in]{geometry}  % Change to margin=0.75in for smaller
```

Then recompile: `latexmk -pdf main.tex`

### Generate without bibliography

To skip bibliography generation:
```bash
pdflatex main.tex
pdflatex main.tex
```

(Omit the `bibtex main.aux` step)

## Output Files

After compilation, your directory will contain:

```
UserManual/
├── main.pdf              ← Final output (read this!)
├── main.tex              ← Master document
├── main.log              ← Compilation log
├── main.aux              ← Auxiliary file
├── main.bbl              ← Bibliography file
├── main.blg              ← Bibliography log
├── main.toc              ← Table of contents
└── ... (other input files)
```

You only need `main.pdf`. The other files are temporary.

To clean up temporary files:
```bash
# macOS/Linux
rm main.aux main.bbl main.blg main.log main.out main.toc

# Or use latexmk
latexmk -c main.tex  # Clean temporary files
latexmk -C main.tex  # Also delete PDF
```

## Reading the PDF

Once `main.pdf` is generated:

### PDF Viewers
- macOS: Preview (built-in) or Adobe Acrobat Reader
- Linux: Okular, Evince, or PDF.js
- Windows: Adobe Acrobat Reader or Windows PDF Viewer

### Navigation Features
1. **Bookmarks panel**: Click bookmark to jump to section
2. **Hyperlinks**: Click blue text to navigate
3. **Search**: Ctrl+F (or Cmd+F on macOS) to search document
4. **Outline**: View Table of Contents in sidebar

### Printing
```bash
# Print to PDF (if you want another copy)
lpr main.pdf

# Or from preview application
```

## Continuous Compilation (Watch Mode)

For development/editing, use latexmk in watch mode:

```bash
latexmk -pdf -pvc main.tex
```

This will:
- Automatically recompile whenever you save `main.tex` or any included file
- Open PDF viewer and refresh it automatically
- Very useful while editing

Press Ctrl+C to stop watching.

## Sharing the Manual

### Share just the PDF
```bash
cp main.pdf ~/Desktop/AxionSR_UserManual.pdf
```

### Share the entire source
```bash
cd ~/Desktop
tar -czf AxionSR_UserManual.tar.gz /Users/samuelwitte/Dropbox/Axion_SR/UserManual/
# Creates AxionSR_UserManual.tar.gz (~150 KB)
```

### Include in repository
```bash
# Copy to code repository
cp -r /Users/samuelwitte/Dropbox/Axion_SR/UserManual ~/AxionSR.jl/docs/

# Or just the PDF
cp main.pdf ~/AxionSR.jl/docs/AxionSR_UserManual.pdf
```

## Expected Compilation Time

- **First compilation**: 30-60 seconds (includes bibliography processing)
- **Subsequent recompilations**: 10-20 seconds
- **Full clean rebuild**: 60 seconds

If compilation takes much longer:
- Check system load: `top` or Activity Monitor
- Close other programs
- Consider older computer may be slower

## Getting Help

If you encounter problems:

1. **Check this guide** for common issues
2. **Check LaTeX log**: Open `main.log` and look for errors
3. **Check file paths**: Ensure all files are in correct directories
4. **Check TeX Live**: Update packages with `tlmgr update --all`
5. **Try clean rebuild**:
   ```bash
   latexmk -C main.tex  # Delete everything
   latexmk -pdf main.tex  # Rebuild from scratch
   ```

## Success Checklist

After compilation:
- [ ] `main.pdf` file exists
- [ ] File size is ~500 KB to 1 MB
- [ ] PDF opens in reader without errors
- [ ] Table of contents is present and clickable
- [ ] Hyperlinks are blue and clickable
- [ ] References in [brackets] appear in bibliography
- [ ] Equations are numbered and can be referenced
- [ ] Page count is 60-80 pages

If all checks pass, you're ready to use the manual!

---

**Last Updated**: November 23, 2025
**Tested On**: macOS 13+ with MacTeX, Ubuntu 22.04 with TeX Live, Windows 10+ with MiKTeX
