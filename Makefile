all: ms/Seedling.bib ms/methods.pdf

ms/methods.pdf: ms/methods.Rmd  ms/Seedling.bib figs/*
	R -e 'system.time(rmarkdown::render("$<", "all"))'

ms/Seedling.bib: ~/Seedling.bib
	cp $< ./ms/

.PHONY: clean
clean:
	rm -f ms/*.tuc \
	ms/*.log
	rm -rf docs/components
