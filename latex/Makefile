#
#******************************************************************************

#
# Locating the root directory
#

fig_dir=figures

out = $(fig_dir)/exp_setup.pdf \
	$(fig_dir)/AFM_loop_z.pdf \
	Andersson_bibs.bib\
	braker_bibs.bib\
	Yufans_bibs.bib



all: afm-cs-main.pdf



$(fig_dir)/exp_setup.pdf: $(fig_dir)/exp_setup.odp
	unoconv --output=tmp.eps -f eps $^
	epstopdf tmp.eps $@
	pdfcrop $@ $@
	rm tmp.eps


$(fig_dir)/AFM_loop_z.pdf: $(fig_dir)/AFM_loop_z.odp
	unoconv --output=tmp.eps -f eps $^
	epstopdf tmp.eps $@
	pdfcrop $@ $@
	rm tmp.eps


afm-cs-main.pdf:afm-cs-main.tex $(out)
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make $<

clean:
	latexmk -CA
