FILE=soptx_cicp

all: 
	pdflatex $(FILE).tex
	bibtex $(FILE).aux
	pdflatex $(FILE).tex
	pdflatex $(FILE).tex
	evince $(FILE).pdf &
	@echo "PDF generated."

clean:
	-rm -fr *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *~ *.tid *.tms *.bak *.blg *.bbl *.gz *.snm *.nav _minted*

view:
	evince $(FILE).pdf &

.PHONY: all clean view
