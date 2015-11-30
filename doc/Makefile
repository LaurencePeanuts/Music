# ----------------------------------------------------------------------------

all: mapping inflation trees

clean:
	-rm -f *.aux *.log *.dvi *.blg *.out *.bbl

# ----------------------------------------------------------------------------

# Paper 1:

mapping: music-mapping.tex macros.tex references.bib
	pdflatex music-mapping.tex
	pdflatex music-mapping.tex
	bibtex music-mapping
	pdflatex music-mapping.tex
	pdflatex music-mapping.tex

# Paper 2:

inflation: music-inflation.tex macros.tex references.bib
	pdflatex music-inflation.tex
	pdflatex music-inflation.tex
	bibtex music-inflation
	pdflatex music-inflation.tex
	pdflatex music-inflation.tex

# Paper 3:

trees: music-trees.tex macros.tex references.bib
	pdflatex music-trees.tex
	pdflatex music-trees.tex
	bibtex music-trees
	pdflatex music-trees.tex
	pdflatex music-trees.tex

# ----------------------------------------------------------------------------
