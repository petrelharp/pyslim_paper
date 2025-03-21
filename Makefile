all : pyslim_chapter.pdf

pyslim_chapter.pdf: pyslim_chapter.tex references.bib

clean:
	rm -f *.pdf *.aux *.bbl *.log

%.pdf : %.tex %.bbl
	while ( pdflatex $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.aux : %.tex
	-pdflatex $<

%.bbl : %.aux
	-bibtex $<

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.pdf : %.md
	pandoc -f markdown -o $@ $<
