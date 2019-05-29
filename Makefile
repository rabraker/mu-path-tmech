
figpath = latex/figures

out = $(figpath)/dna_og.svg


all:block_diagrams $(out)

block_diagrams:
	cd latex && $(MAKE)


$(figpath)/dna_og.svg:plot_bptv_vs_bp.m
	./matlab_wrapper.py $<
