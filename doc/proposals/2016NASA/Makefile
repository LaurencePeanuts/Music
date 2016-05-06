# ----------------------------------------------------------------------------

all: proposal

clean:
	-rm -f *.aux *.log *.dvi *.blg *.out *.bbl

# ----------------------------------------------------------------------------

proposal: prodes.tex prodes.bib
	pdflatex prodes.tex
	pdflatex prodes.tex
	bibtex prodes
	pdflatex prodes.tex
	pdflatex prodes.tex

# ----------------------------------------------------------------------------
