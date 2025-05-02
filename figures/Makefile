.PHONY : all

INKSCAPE := $(shell command -v inkscape 2> /dev/null)

# NOTE: all layers must be visible in pedigrees.ink.svg for this to work.

all: pedigree0.pdf pedigree1.pdf pedigree2.pdf pedigree_recapitate.pdf pedigree_simplify.pdf pedigree_mutate.pdf

pedigree0.svg : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		ts_indivs simp_indivs extra_indivs extra_pedigree pedigree simp_pedigree \
		remember_pedigree chromosomes simplify_removes extra_chromosomes label0 \
		> $@

pedigree1.svg : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		time_label ts_indivs simp_indivs pedigree simp_pedigree \
		chromosomes simplify_removes label1 \
		> $@

pedigree2.svg : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		time_label ts_indivs simp_indivs pedigree simp_pedigree \
		chromosomes simplify_removes sample_labels node_labels \
		> $@

pedigree_recapitate.svg : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		tskit_time ts_indivs simp_indivs pedigree simp_pedigree recap_pedigree \
		chromosomes recap_chromosomes simplify_removes recap_label \
		> $@

pedigree_simplify.svg : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		tskit_time ts_indivs pedigree recap_pedigree \
		chromosomes recap_chromosomes \
		> $@

pedigree_mutate.svg : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		ts_indivs pedigree chromosomes \
		tskit_time ts_indivs pedigree recap_pedigree \
		chromosomes recap_chromosomes mutations \
		> $@

pedigree_remember.svg : pedigrees.ink.svg
	./export-layers-svg.sh $< \
		ts_indivs simp_indivs remember_indivs pedigree simp_pedigree remember_pedigree \
		chromosomes simplify_removes sample_labels remember_labels \
		> $@

## 
# Graphics whatnot

%.png : %.pdf
ifndef INKSCAPE
	echo "inkscape not available: not regenerating figures"
else
	inkscape $<  --export-type="png" --export-dpi=300
endif

%.png : %.svg
ifndef INKSCAPE
	echo "inkscape not available: not regenerating figures"
else
	inkscape $<  --export-type="png" --export-dpi=300
endif

%.pdf : %.svg
ifndef INKSCAPE
	echo "inkscape not available: not regenerating figures"
else
	inkscape $<  --export-type="pdf"
endif

clean :
	-rm pedigree*.png $(filter-out pedigrees.ink.svg,$(wildcard *.svg))
