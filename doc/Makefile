main.pdf:	library2.bib CodeStructure.tex Equations3.tex ExistingModels.tex ExistingModels4.tex Introduction.tex main.tex MatlabCode.tex ModelDevelopment.tex Results2.tex Results3.tex ToDo2.tex
	pdflatex main
	bibtex main
	pdflatex main
	makeindex -s main.ist -t main.glg -o main.gls main.glo
	makeindex -s main.ist -t main.slg -o main.syi main.syg
	pdflatex main

	
