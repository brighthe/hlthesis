
all: phd_research_proposal.pdf \

subspace_correct.pdf: subspace_correct.tex
	xelatex -shell-escape subspace_correct.tex
	#evince subspace_correct.pdf&

phd_research_proposal.pdf: phd_research_proposal.tex
	xelatex -shell-escape phd_research_proposal.tex
	bibtex phd_research_proposal.aux
	xelatex -shell-escape phd_research_proposal.tex
	xelatex -shell-escape phd_research_proposal.tex
	evince phd_research_proposal.pdf&

gravity_wave.pdf: gravity_wave.tex
	xelatex -shell-escape gravity_wave.tex
	xelatex -shell-escape gravity_wave.tex
	evince gravity_wave.pdf&

infsup.pdf: infsup.tex
	xelatex -shell-escape infsup.tex
	xelatex -shell-escape infsup.tex
	#evince infsup.pdf&

hodge_wave_equation.pdf: hodge_wave_equation.tex
	xelatex -shell-escape hodge_wave_equation.tex
	xelatex -shell-escape hodge_wave_equation.tex
	#evince hodge_wave_equation.pdf&

.PHONY:clean  
clean:
	-rm -fr *.xml *.bcf *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *~ *.tid *.tms *.pdf *.bak *.blg *.bbl *.gz *.snm *.nav _minted*

