######################################################################
##########                                                  ##########
##########        Start Date: April.20 2019                 ##########
##########        Author: Wei Si                            ##########
##########        Email: siwei@smail.xtu.edu.cn      		##########
##########        Modified for brightPhD                    ##########
##########                                                  ##########
######################################################################

#### Compile Latex file;
### all: file   ('all' depends on 'file');
### .PHONY: object  ('object' is pointed a phony object);
### @  (display the code in the terminal);
### -  (ignore the error of the command);

#### If you use this Makefile, please replace 'file' by your file name;

file = brightPhD


.PHONY: all
all:
ifeq ($(file).tex, $(wildcard $(file).tex))
	@-xelatex $(file).tex
	@-bibtex $(file).aux
	@-xelatex $(file).tex
	@-xelatex $(file).tex
	clear
else
	clear
	@echo "WARNING: Lack tex files"
endif


.PHONY: vim
vim:
ifeq ($(file).tex, $(wildcard $(file).tex))
	clear
	@vim $(file).tex
else
	clear
	@echo "WARNING: Lack tex files"
endif


.PHONY: view
view:
ifeq ($(file).pdf, $(wildcard $(file).pdf))
	clear
	@evince $(file).pdf &
else
	clear
	@echo "WARNING: Lack pdf files"
endif


.PHONY: clean
clean:
	clear
	@-rm -f *.ps *.dvi *.aux *.toc *.idx *.ind *.ilg *.log *.out *~ *.tid *.tms *.bak *.blg *.bbl *.nav *.snm *.vrb *.thm *.spl *.soc *.bcf *.xml *.synctex.gz
	@-rm -f body/*.aux
	@-rm -f body/*/*.aux
	@-rm -f reference/*.aux
	@-rm -f preface/*.aux
	@echo "Hint: Clean Out"


.PHONY: cleanall
cleanall: clean
	@-rm -f $(file).pdf
	@echo "Hint: Clean All (including PDF)"


.PHONY: help
help:
	clear
	@echo "==================== A Common Makefile for Latex ===================="
	@echo " Copyright (Latex)	2019	Wei Si"
	@echo " Modified for brightPhD project"
	@echo " The following targets are support"
	@echo
	@echo " all				- (== make) compile the document"
	@echo " vim				- edit the main file"
	@echo " view				- view the pdf file"
	@echo " clean			- clean temporary files from compiling"
	@echo " cleanall			- clean all files including PDF"
	@echo "========================== Version 1.1 =========================="