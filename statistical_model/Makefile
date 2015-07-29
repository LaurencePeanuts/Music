# ----------------------------------------------------------------------------

all: ctu2015

clean:
	-rm -f *.aux *.log *.dvi *.blg *.out *.bbl ctu2015.pdf

# ----------------------------------------------------------------------------

# Notes from CtU2015:

ctu2015: ctu2015.tex macros.tex references.bib
	pdflatex ctu2015.tex
	pdflatex ctu2015.tex
	bibtex ctu2015
	pdflatex ctu2015.tex
	pdflatex ctu2015.tex

# ----------------------------------------------------------------------------
