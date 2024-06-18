
all: Timo.pdf \

Timo.pdf: Timo.tex
	xelatex -shell-escape Timo.tex
	xelatex -shell-escape Timo.tex
	evince Timo.pdf&

.PHONY:clean  
clean:
	-rm -fr  */*.vrb */*.bcf */*.mtc */*.maf */*.mtc* */*.xml */*.ps */*.dvi */*.aux */*.toc */*.idx */*.ind */*.ilg */*.log */*.out */*~ */*.tid */*.tms */*.bak */*.blg */*.bbl */*.gz */*.snm */*.nav */_minted*

