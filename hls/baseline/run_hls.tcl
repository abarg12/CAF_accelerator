# Vitis HLS script for the CAF baseline.
#
#   vitis-run --mode hls --tcl run_hls.tcl
#
# Runs C simulation, C synthesis, and C/RTL co-simulation in caf_baseline.proj.
#
# The design is built with 8 Doppler bins instead of the full 256, because
# co-simulating all 256 bins takes hours in the RTL simulator. Every Doppler
# bin does the same work, so the full design's latency is about 256 / 8 = 32
# times the measured one. Set CAF_FREQ_BINS to build a different size, e.g.
# CAF_FREQ_BINS=256 for the full design's resources (expect a long cosim).

set freq_bins 8
if {[info exists ::env(CAF_FREQ_BINS)]} {
	set freq_bins $::env(CAF_FREQ_BINS)
}
set cflags "-DNUM_OF_FREQ_BINS=$freq_bins"

open_project caf_baseline.proj -reset
add_files caf.cpp -cflags $cflags
add_files fft.cpp -cflags $cflags
add_files -tb caf_test.cpp -cflags $cflags
set_top cross_ambiguity

open_solution solution -reset

# PYNQ-Z2 (Zynq-7020), 100 MHz clock
set_part xc7z020clg400-1
create_clock -period 10

csim_design
csynth_design
cosim_design -O

exit
