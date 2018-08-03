name = cslvr_book

all:
	sh ./build.sh $(name)

final:
	sh ./build.sh $(name) --final

chapter:
	BIBINPUTS="src:" 
	sh ./build.sh $(name) --chapter $(chapter)

clean: 
	rm texput.log $(name).aux $(name).bbl $(name).blg $(name).idx $(name).ilg $(name).ind $(name).loa $(name).lof $(name).log $(name).lot $(name).out $(name).pdf $(name).run.xml $(name).toc $(name)-blx.bib



