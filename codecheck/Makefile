all: codecheck.pdf

codecheck.pdf: codecheck.Rmd codecheck-preamble.sty ../codecheck.yml
	Rscript -e 'rmarkdown::render("$<", "pdf_document")'

.PHONY: run clean

clean:
	rm -fr codecheck.pdf outputs


