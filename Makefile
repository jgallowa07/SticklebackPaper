.PHONY: all, clean

# don't delete these
.PRECIOUS: Stickleback_Paper.aux Stickleback_Paper.bbl

all: Stickleback_Paper.pdf

submission : Stickleback_Paper_text.pdf submission/resubmission_cover.pdf Stickleback_Paper_responses.pdf

Stickleback_Paper_text.pdf : Stickleback_Paper.pdf
	pdfjam --outfile $@ $< 1-27

submission/resubmission_cover.pdf : Stickleback_Paper.pdf
	pdfjam --outfile $@ $< 28

Stickleback_Paper_responses.pdf : Stickleback_Paper.pdf
	pdfjam --outfile $@ $< 29-

clean: 
	-rm *.aux *.log *.lof *.lot *.fff *.ttt *.out *.bbl *.blg

%.pdf : %.tex %.bbl
	while ( pdflatex $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.aux : %.tex
	-pdflatex $<

%.bbl : %.aux
	bibtex $<

%.html : %.md
	Rscript -e "templater::render_template(md.file='$<', output='$@')"

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

%.eps : %.pdf
	inkscape --without-gui --export-eps=$@ $<
