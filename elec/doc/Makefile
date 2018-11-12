# Thesis make file

LATEX = pdflatex

BIBTEX = bibtex

all: cardioid_manual.pdf

cardioid_manual.pdf: Cardioid.tex
	$(LATEX) Cardioid
	$(LATEX) Cardioid
	$(LATEX) Cardioid

clean:
	rm -f *.dvi
	rm -f *.bbl
	rm -f *.toc
	rm -f *.blg
	rm -f *.aux
	rm -f *.ps
	rm -f *.pdf
