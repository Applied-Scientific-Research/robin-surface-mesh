# type "make" command in Unix to create asme2e.pdf file 
#ROOT=robin
ROOT=preprint

all:
	lualatex $(ROOT)
	bibtex $(ROOT)
	lualatex $(ROOT)
	#lualatex $(ROOT)

clean:
	(rm -rf *.ps *.log *.dvi *.aux *.*% *.lof *.lop *.lot *.toc *.idx *.ilg *.ind *.bbl *blg)
