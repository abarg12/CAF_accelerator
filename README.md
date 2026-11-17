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
| `baseline/signal_generation.py` | Simulated FM (or white noise) transmit signal, a delayed, Doppler shifted, noisy receive signal, and stationary clutter |
| `baseline/run_baseline.py` | Plots a range-Doppler map and times the CAF from 8x8 to 256x256 bins |
| `tests/test_caf_baseline.py` | Checks that the CAF recovers known targets and clutter, and that the FM signal is band limited |

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

### Simulated scenario

`run_baseline.py` builds a scenario like the paper's (Sections 2.5 and 2.6,
Figure 11):

- **Transmit signal:** a simulated FM broadcast signal (75 kHz deviation,
  53 kHz stereo multiplex baseband, about 256 kHz wide). The paper records a
  real FM station with an SDR instead.
- **Target:** one moving target at 10 km with a 12 kHz Doppler shift.
- **Stationary clutter:** four zero Doppler returns at 3, 16.5, 24, and
  31.5 km.
- **Noise:** complex Gaussian noise at 0 dB SNR relative to the target. The
  256 sample CAF has about 24 dB of processing gain, so the target stays well
  above the noise.

The sample rate is 1 MHz, so each range bin is 150 m and 256 bins cover
38 km. The Doppler extent is 64 kHz (+/-32 kHz, 250 Hz per bin), set apart
from `Fs/N` as in Eq. 1. The true Doppler resolution is `Fs/N` = 3.9 kHz, so
a target spreads over several Doppler bins.

Two plots are saved. `range_doppler_map.png` shows `|CAF|` in dB relative
to the peak, drawn as one flat square per bin with no interpolation.
`range_doppler_map_cuts.png` shows the row (Doppler cut) and the column (range
cut) of the map through the detected target, drawn as steps, one per bin.
The map looks smooth because neighboring bins hold similar values: a target
spans about 16 Doppler bins (3.9 kHz resolution / 250 Hz bins) and about 4
range bins (the FM signal's ~590 m range resolution / 150 m bins). The target is detected as the
strongest return outside the zero Doppler clutter band. The FPGA core and
`calculate_cross_ambiguity` return the real part, as in the paper. The
magnitude is used only for plotting, because the real part swings with the
phase of each return.

With only 256 samples the FM signal's own sidelobes sit roughly 10 to 15 dB
below a peak, which is why the background is not flat. Clutter much stronger
than the target would hide it in those sidelobes. Real passive radars remove
clutter before the CAF to deal with this.

The tests keep using the white noise transmit signal, because its one bin
peak makes exact delay and Doppler checks easy.

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
| `hls/run_hls.tcl` | C simulation, C synthesis, and C/RTL co-simulation for the PYNQ-Z2 (xc7z020clg400-1) at 100 MHz |

### Running it

From the `hls/` folder:

    vitis-run --mode hls --tcl run_hls.tcl

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
