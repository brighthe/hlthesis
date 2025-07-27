
all: cbtxsPHD.pdf \

cbtxsPHD.pdf: cbtxsPHD.tex
	xelatex -shell-escape cbtxsPHD.tex
	bibtex cbtxsPHD
	xelatex -shell-escape cbtxsPHD.tex
	xelatex -shell-escape cbtxsPHD.tex
	#evince cbtxsPHD.pdf&


.PHONY:clean  
clean:
	-rm -fr *.aux *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *~ *.tid *.tms *.bak *.blg *.bbl *.gz *.snm *.nav _minted*

