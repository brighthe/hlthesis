
all: main.pdf \

main.pdf: main.tex
	pdflatex -shell-escape main.tex
	#bibtex main.aux
	pdflatex -shell-escape main.tex
	pdflatex -shell-escape main.tex
	#evince main.pdf&

test.pdf: test.tex
	xelatex -shell-escape test.tex
	#bibtex test.aux
	xelatex -shell-escape test.tex
	xelatex -shell-escape test.tex
	#evince test.pdf&

mid_term_report.pdf: mid_term_report.tex
	xelatex -shell-escape mid_term_report.tex
	#bibtex mid_term_report.aux
	xelatex -shell-escape mid_term_report.tex
	xelatex -shell-escape mid_term_report.tex
	#evince mid_term_report.pdf&

stabfree.pdf: stabfree.tex
	xelatex -shell-escape stabfree.tex
	#bibtex stabfree.aux
	xelatex -shell-escape stabfree.tex
	xelatex -shell-escape stabfree.tex
	#evince stabfree.pdf&

defense.pdf: defense.tex
	xelatex -shell-escape defense.tex
	#xelatex -shell-escape defense.tex
	evince defense.pdf&


.PHONY:clean  
clean:
	-rm -fr *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *~ *.tid *.tms *.pdf *.bak *.blg *.bbl *.gz *.snm *.nav _minted*

