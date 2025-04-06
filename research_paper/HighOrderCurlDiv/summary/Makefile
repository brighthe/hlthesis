FILE=summary

all: 
	pdflatex $(FILE).tex
	pdflatex $(FILE).tex
	evince $(FILE).pdf &
	@echo "PDF generated."

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb

view:
	evince $(FILE).pdf &

.PHONY: all clean view
