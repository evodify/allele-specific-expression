MAIN=tutorial

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex
	latexmk -pdf $(MAIN).tex

clean:
	rm -f $(MAIN).aux $(MAIN).bbl $(MAIN).blg $(MAIN).dvi $(MAIN).fdb_latexmk $(MAIN).fls $(MAIN).log $(MAIN).out
