# type "make" command in Unix to create asme2e.pdf file 
ROOT=robin

all:
	lualatex $(ROOT)
	bibtex $(ROOT)
	lualatex $(ROOT)
	#lualatex $(ROOT)

abstract:
	latex abstract
	dvips -o abstract.ps abstract
	ps2pdf abstract.ps abstract.pdf

clean:
	(rm -rf *.ps *.log *.dvi *.aux *.*% *.lof *.lop *.lot *.toc *.idx *.ilg *.ind *.bbl *blg)
