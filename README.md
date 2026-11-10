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
