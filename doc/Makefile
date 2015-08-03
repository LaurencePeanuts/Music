# ----------------------------------------------------------------------------

all: music

clean:
	-rm -f *.aux *.log *.dvi *.blg *.out *.bbl

# ----------------------------------------------------------------------------

# Notes from CtU2015:

ctu2015: ctu2015.tex macros.tex references.bib
	pdflatex ctu2015.tex
	pdflatex ctu2015.tex
	bibtex ctu2015
	pdflatex ctu2015.tex
	pdflatex ctu2015.tex

# Working report:

music: music.tex references.bib
	pdflatex music.tex
	pdflatex music.tex
	bibtex ctu2015
	pdflatex music.tex
	pdflatex music.tex

# ----------------------------------------------------------------------------
