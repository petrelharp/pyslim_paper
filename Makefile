.PHONY : figures

all : pyslim_chapter.pdf

pyslim_chapter.pdf : pyslim_chapter.tex references.bib figures

figures :
	$(MAKE) -C figures

clean:
	rm -rf *.pdf *.aux *.bbl *.log _minted

%.pdf : %.tex %.bbl
	while ( pdflatex --shell-escape $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.aux : %.tex
	-pdflatex --shell-escape $<

%.bbl : %.aux
	-bibtex $<

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.pdf : %.md
	pandoc -f markdown -o $@ $<
