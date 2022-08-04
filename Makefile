GIT = c31e6ea2
METHODS = ms/methods
all: $(METHODS).pdf $(METHODS).docx ms/methods-diff$(GIT).pdf
diff: ms/methods-diff$(GIT).pdf

$(METHODS).pdf: $(METHODS).qmd
	quarto render $< --to pdf

$(METHODS).docx: $(METHODS).qmd
	quarto render $< --to docx

ms/methods-diff$(GIT).pdf: ms/methods.tex
	latexdiff-vc --git --flatten --force -r $(GIT) $^ ; \
	cd ms; pdflatex methods-diff$(GIT).tex

.PHONY: clean
clean:
	rm -f ms/*.tuc \
	ms/*.log
	rm -rf docs/components
