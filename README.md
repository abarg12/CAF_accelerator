# CAF Accelerator

FPGA acceleration of the radar Cross Ambiguity Function (CAF) with Vitis HLS,
compared against a CPU baseline.

## CPU baseline

The Python baseline in `baseline/` is the software reference model. It follows
the prototype in Figure 4 of Long et al., "Hardware-in-the-Loop Simulation with
PYNQ for Radar Applications" (GNU Radio Conference 2025), and computes

    X(v) = IFFT( conj(FFT(x)) * FFT(y * exp(2j*pi*v*m / Fs)) ) / N

for every Doppler bin `v`, where `x` is the transmit signal and `y` is the
receive signal.

| File | Purpose |
| --- | --- |
| `baseline/caf_baseline.py` | The CAF itself: frequency shifts, FFT cross correlation, range-Doppler map |
| `baseline/signal_generation.py` | Simulated transmit signal and a delayed, Doppler shifted, noisy receive signal |
| `baseline/run_baseline.py` | Plots a range-Doppler map and times the CAF from 8x8 to 256x256 bins |
| `tests/test_caf_baseline.py` | Checks that the CAF recovers a known target delay and Doppler |

### Running it

The Python environment is managed with `venv`. Create it once from the repo root:

    python3 -m venv .venv
    source .venv/bin/activate        # Windows: .venv\Scripts\activate
    pip install -r requirements.txt

Then, with the environment activated:

    python baseline/run_baseline.py
    python -m pytest tests

Run `deactivate` to leave the environment. The `.venv/` folder is ignored by
git, so each machine builds its own from `requirements.txt`.

The range-Doppler map returned by `calculate_cross_ambiguity` has shape
`(num_of_samples, num_of_freq_bins)`: rows are delay (range) bins and columns
are Doppler bins. `delay_axis` and `doppler_axis` give the labels for each.

## HLS baseline

The Vitis HLS baseline in `hls/` is a direct C++ port of the Python reference
model, with no HLS pragmas. It gives the starting latency and resource numbers
that later optimizations are measured against.

| File | Purpose |
| --- | --- |
| `hls/caf.cpp` | Top function `cross_ambiguity`: frequency shift, FFT cross correlation, range-Doppler map |
| `hls/fft.cpp` | Radix-2 FFT and inverse FFT |
| `hls/caf.h` | Sizes (256 samples, 256 Doppler bins), data type, and prototypes |
| `hls/caf_test.cpp` | Self checking testbench |
| `hls/run_hls.tcl` | C simulation and C synthesis for the PYNQ-Z2 (xc7z020clg400-1) at 100 MHz |

### Running it

With Vitis HLS, from the `hls/` folder:

    vitis_hls -f run_hls.tcl                      # Vitis HLS 2023.2 and earlier
    vitis-run --mode hls --tcl run_hls.tcl        # Vitis 2024.1 and later

The testbench also runs without Vitis:

    g++ -o caf_test caf.cpp fft.cpp caf_test.cpp && ./caf_test

## References

- C. D. Long, T. Brant, L. Brand, and J. D. Ouellette, "Hardware-in-the-Loop
  Simulation with PYNQ for Radar Applications," Proceedings of the 15th GNU
  Radio Conference, 2025. The Python baseline and the CAF structure follow
  this paper.
- R. Kastner, J. Matai, and S. Neuendorffer, "Parallel Programming for FPGAs,"
  https://github.com/KastnerRG/pp4fpgas. The FFT in `hls/fft.cpp` is adapted
  from the `fft_sw` example in this repository, which is licensed under
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The changes are
  listed at the top of that file.
