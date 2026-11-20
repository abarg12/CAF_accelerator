# Vitis HLS script for the CAF baseline.
#
#   vitis-run --mode hls --tcl run_hls.tcl
#
# 1. Runs C simulation and C synthesis of the full design (256 Doppler bins)
#    in caf_baseline.proj.
# 2. Runs C synthesis and C/RTL co-simulation of a smaller design with 8
#    Doppler bins in caf_cosim_8bins.proj. Co-simulating all 256 bins takes
#    hours in the RTL simulator. Every Doppler bin does the same work, so the
#    full design's latency is about 256 / 8 = 32 times the measured one.
#
# Set CAF_COSIM_BINS to co-simulate a different number of bins (even), or to
# 0 to skip co-simulation:
#
#   CAF_COSIM_BINS=4 vitis-run --mode hls --tcl run_hls.tcl

set cosim_bins 8
if {[info exists ::env(CAF_COSIM_BINS)]} {
	set cosim_bins $::env(CAF_COSIM_BINS)
}

# Create a project for the PYNQ-Z2 (Zynq-7020) at 100 MHz. cflags is passed
# to the design and the testbench, e.g. to change NUM_OF_FREQ_BINS.
proc open_caf_project {project cflags} {
	open_project $project -reset
	add_files caf.cpp -cflags $cflags
	add_files fft.cpp -cflags $cflags
	add_files -tb caf_test.cpp -cflags $cflags
	set_top cross_ambiguity

	open_solution solution -reset
	set_part xc7z020clg400-1
	create_clock -period 10
}

# Full design
open_caf_project caf_baseline.proj ""
csim_design
csynth_design
close_project

# Smaller design for co-simulation
if {$cosim_bins > 0} {
	open_caf_project caf_cosim_${cosim_bins}bins.proj "-DNUM_OF_FREQ_BINS=$cosim_bins"
	csynth_design
	cosim_design -O
	close_project
}

exit
