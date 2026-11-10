# Vitis HLS script for the CAF baseline.
#
#   vitis_hls -f run_hls.tcl                      (Vitis HLS 2023.2 and earlier)
#   vitis-run --mode hls --tcl run_hls.tcl        (Vitis 2024.1 and later)

open_project caf_baseline.proj -reset
add_files caf.cpp
add_files fft.cpp
add_files -tb caf_test.cpp

set_top cross_ambiguity

open_solution solution -reset

# PYNQ-Z2 (Zynq-7020), 100 MHz clock
set_part xc7z020clg400-1
create_clock -period 10

csim_design
csynth_design

exit
