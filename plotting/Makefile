all : merged.muts.svg

out.trees : back-and-forth.slim
	slim -s 123 $<

%.muts.svg : %.tsz
	python back-and-forth.py $< 
