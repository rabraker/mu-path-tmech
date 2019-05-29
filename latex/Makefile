#
#******************************************************************************

#
# Locating the root directory
#

fig_dir=figures

out = $(fig_dir)/exp_setup.pdf \
	$(fig_dir)/AFM_loop_z.pdf



all: $(out)


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
