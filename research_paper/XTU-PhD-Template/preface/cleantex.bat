@echo off
echo clean started at (%TIME%)
echo written by Xiaoqiang Yue
del ctextemp*
del *.aux
del *.log
del *.nav
del *.out
del *.djs
del *.out.bak
del *.snm
del *.ps
del *.dvi
del *.synctex.gz*
del *.idx
del *.toc
del *.tex.bak
del *.blg
del *.bbl
del *.ilg
del *.ind
del *.tex.sav
del *.sty.bak
del *.sty.sav
del *.fls
del *.spl
del *.fdb_latexmk
echo clean finished at (%TIME%)