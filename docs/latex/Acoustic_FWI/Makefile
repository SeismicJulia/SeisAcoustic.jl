textpdf :
	pdflatex paper
	bibtex paper
	pdflatex paper
	pdflatex paper
	make tidy

textps :
	latex paper
	bibtex paper
	latex paper
	latex paper
	dvips -t letter -P pdf paper -o
	ps2pdf paper.ps

clean : tidy
	rm -f paper.pdf

tidy :
	rm -f *.aux
	rm -f *~
	rm -f *.bbl
	rm -f *.blg
	rm -f *.log
	rm -f *.ps
	rm -f *.dvi
	rm -f *.bak
	rm -f *.fff
	rm -f *.lof


