
figpath = latex/figures

out = $(figpath)/dna_og.svg \
	$(figpath)/MIMO_CL_uxuy.svg \
	$(figpath)/z_control_design.svg \
	$(figpath)/with_wo_ff.svg       \
	$(figpath)/baseline_errors_aligned_1Hz.svg \
	$(figpath)/colorbar.svg

# Stupid overleaf always changes the file permissions of matlab_wrapper.py to
# 644 (non-executable), so we cant run it as ./matlab_wrapper.py
#

all: $(out)
	cd latex && $(MAKE)


$(figpath)/dna_og.svg:plot_bptv_vs_bp.m functions/Core_Nesterov_mine.m functions/NESTA_mine.m
	python matlab_wrapper.py $<

$(figpath)/with_wo_ff.svg:plot_with_without_FF_loop_cs20ngflat.m
	python matlab_wrapper.py $<

$(figpath)/MIMO_CL_uxuy.svg:build_controllers_with_xyFF_DxDy.m
	python matlab_wrapper.py $<


$(figpath)/z_control_design.svg:plot_z_axis_control_design.m
	python matlab_wrapper.py $<

$(figpath)/baseline_errors_aligned_1Hz.svg:plot_baseline_raster_analysis_1Hz.m
	python matlab_wrapper.py $<

$(figpath)/colorbar.svg:make_colorbar.m
	python matlab_wrapper.py $<
