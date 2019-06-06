
figpath = latex/figures

out = $(figpath)/dna_og.svg \
	$(figpath)/MIMO_CL_uxuy.svg \
	$(figpath)/z_control_design.svg


all:block_diagrams $(out)

block_diagrams:
	cd latex && $(MAKE)


$(figpath)/dna_og.svg:plot_bptv_vs_bp.m functions/Core_Nesterov_mine.m functions/NESTA_mine.m
	./matlab_wrapper.py $<


$(figpath)/MIMO_CL_uxuy.svg:build_controllers_with_xyFF_DxDy.m
	./matlab_wrapper.py $<


$(figpath)/z_control_design.svg:plot_z_axis_control_design.m
	./matlab_wrapper.py $<
