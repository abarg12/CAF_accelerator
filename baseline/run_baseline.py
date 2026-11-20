"""
Run the CAF CPU baseline: plot a range-Doppler map and time the computation.

    python baseline/run_baseline.py              # plot + benchmark
    python baseline/run_baseline.py --no-plot    # benchmark only

The simulated scenario follows Sections 2.5 and 2.6 and Figure 11 of the paper:
an FM broadcast signal as the transmit signal, one moving target near 10 km,
stationary clutter at several ranges, and additive noise.

The benchmark sweeps square range-Doppler maps from 8x8 to 256x256, the same
sizes reported in the paper, so the numbers can be compared against the FPGA
core later. The run time does not depend on the signal content.
"""

import argparse
import os
import time

import numpy as np

from caf_baseline import cross_ambiguity

SPEED_OF_LIGHT = 299_792_458.0  # m/s
SAMPLE_RATE = 1_000_000  # Hz

# Total Doppler extent R from Eq. 1, centered on 0 Hz. The Doppler step is
# R / num_of_freq_bins. This is finer than the Doppler resolution of
# sample_rate / num_of_samples (3.9 kHz for 256 samples), so a target spreads
# over several Doppler bins, as in Figure 11 of the paper.
DOPPLER_EXTENT_HZ = 64_000

TARGET_RANGE_METERS = 10_000
TARGET_DOPPLER_HZ = 12_000
SIGNAL_TO_NOISE_RATIO_DB = 0  # target return power relative to noise power

# Stationary clutter: (range in meters, amplitude relative to the target)
CLUTTER = [(3_000, 0.6), (16_500, 0.5), (24_000, 0.55), (31_500, 0.4)]

BENCHMARK_BIN_SIZES = [8, 16, 32, 64, 128, 256]
BENCHMARK_REPEATS = 20


def make_signals(num_of_samples, random_generator):
    """Simulated transmit and receive signals with one moving target and clutter.

    Returns (transmit_signal, receive_signal, target_delay, clutter_delays),
    with delays in samples. Delays are circular, to line up with the circular
    correlation of the FFT based CAF, and wrap around for maps too small to
    hold the full range.
    """
    sample_times = np.arange(num_of_samples) / SAMPLE_RATE

    def range_to_delay(range_meters):
        return round(2 * range_meters / SPEED_OF_LIGHT * SAMPLE_RATE) % num_of_samples

    # Transmit: an FM broadcast signal at complex baseband. The message stands
    # in for the FM stereo multiplex baseband (audio, pilot, and stereo
    # subcarrier out to about 53 kHz): white noise low pass filtered to 53 kHz
    # and scaled to a peak of 1. It is frequency modulated with the broadcast
    # deviation of 75 kHz, so by Carson's rule the signal is about
    # 2 * (75 + 53) = 256 kHz wide.
    noise_spectrum = np.fft.fft(random_generator.standard_normal(num_of_samples))
    frequencies = np.fft.fftfreq(num_of_samples, d=1 / SAMPLE_RATE)
    noise_spectrum[np.abs(frequencies) > 53_000] = 0
    message = np.fft.ifft(noise_spectrum).real
    message /= np.max(np.abs(message))
    transmit_signal = np.exp(2j * np.pi * 75_000 * np.cumsum(message) / SAMPLE_RATE)

    # Receive: the target return (delayed and Doppler shifted), plus noise,
    # plus zero Doppler returns from stationary clutter
    target_delay = range_to_delay(TARGET_RANGE_METERS)
    receive_signal = np.roll(transmit_signal, target_delay) * np.exp(
        2j * np.pi * TARGET_DOPPLER_HZ * sample_times
    )
    noise_power = 10 ** (-SIGNAL_TO_NOISE_RATIO_DB / 10)
    receive_signal += np.sqrt(noise_power / 2) * (
        random_generator.standard_normal(num_of_samples)
        + 1j * random_generator.standard_normal(num_of_samples)
    )
    clutter_delays = [range_to_delay(r) for r, _ in CLUTTER]
    for delay, (_, amplitude) in zip(clutter_delays, CLUTTER):
        receive_signal += amplitude * np.roll(transmit_signal, delay)

    return (
        transmit_signal.astype(np.complex64),
        receive_signal.astype(np.complex64),
        target_delay,
        clutter_delays,
    )


def run_benchmark(random_generator):
    print(f"{'bins':>10} {'mean (ms)':>12} {'min (ms)':>12}")
    for bin_size in BENCHMARK_BIN_SIZES:
        transmit_signal, receive_signal, _, _ = make_signals(bin_size, random_generator)
        freq_step = DOPPLER_EXTENT_HZ / bin_size
        run_times_ms = []
        for _ in range(BENCHMARK_REPEATS):
            start_time = time.perf_counter()
            cross_ambiguity(transmit_signal, receive_signal, bin_size, freq_step, SAMPLE_RATE).real
            run_times_ms.append((time.perf_counter() - start_time) * 1e3)
        print(
            f"{f'{bin_size}x{bin_size}':>10} "
            f"{np.mean(run_times_ms):>12.3f} {np.min(run_times_ms):>12.3f}"
        )


def plot_range_doppler_map(random_generator, output_path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    num_of_samples = 256
    num_of_freq_bins = 256
    freq_step = DOPPLER_EXTENT_HZ / num_of_freq_bins
    dynamic_range_db = 20

    transmit_signal, receive_signal, target_delay, clutter_delays = make_signals(
        num_of_samples, random_generator
    )
    caf = cross_ambiguity(transmit_signal, receive_signal, num_of_freq_bins, freq_step, SAMPLE_RATE)

    # Plot |CAF| in dB relative to the peak. The FPGA core outputs the real
    # part, which swings with the phase of each return; the magnitude does
    # not, so it is easier to read.
    magnitude = np.abs(caf)
    map_db = 20 * np.log10(magnitude / magnitude.max() + 1e-12)
    doppler_khz = (np.arange(num_of_freq_bins) - num_of_freq_bins // 2 + 1) * freq_step / 1e3
    range_km = SPEED_OF_LIGHT * np.arange(num_of_samples) / SAMPLE_RATE / 2 / 1e3

    # Detect the moving target as the strongest return outside the zero
    # Doppler clutter band (one Doppler resolution cell on each side of 0 Hz).
    clutter_band = np.abs(doppler_khz) < SAMPLE_RATE / num_of_samples / 1e3
    moving_map_db = np.where(clutter_band, -np.inf, map_db)
    peak_row, peak_column = np.unravel_index(np.argmax(moving_map_db), map_db.shape)

    # Figure 1: the range-Doppler map. pcolormesh draws one flat square per
    # bin with no interpolation. The figure is large enough that each of the
    # 256 x 256 bins is several pixels wide.
    map_figure, map_axis = plt.subplots(figsize=(9, 8))
    image = map_axis.pcolormesh(
        doppler_khz, range_km, map_db, vmin=-dynamic_range_db, vmax=0, shading="nearest",
        cmap="viridis", antialiased=False,
    )
    map_axis.annotate(
        "Target", (doppler_khz[peak_column], range_km[peak_row]),
        xytext=(14, 16), textcoords="offset points", color="white", fontsize=9,
        arrowprops={"arrowstyle": "->", "color": "white", "lw": 0.8},
    )
    map_axis.scatter(
        np.zeros(len(clutter_delays)), range_km[clutter_delays], s=60, facecolors="none",
        edgecolors="white", linewidths=0.9, label="Stationary clutter (placed)",
    )
    map_axis.legend(loc="upper right", fontsize=8, framealpha=0.6)
    map_axis.set_xlabel(f"Doppler shift (kHz), {freq_step:g} Hz per bin")
    map_axis.set_ylabel(f"Range (km), {range_km[1] * 1e3:.0f} m per bin")
    map_axis.set_title("CAF CPU baseline: range-Doppler map")
    map_figure.colorbar(image, ax=map_axis, label="|CAF| (dB relative to peak)")
    map_figure.savefig(output_path, dpi=200, bbox_inches="tight")

    # Figure 2: one row and one column of the map through the detected target,
    # drawn as steps so that each flat segment is one bin.
    stem, extension = os.path.splitext(output_path)
    cuts_path = f"{stem}_cuts{extension or '.png'}"
    cuts_figure, (doppler_cut_axis, range_cut_axis) = plt.subplots(2, 1, figsize=(9, 7))

    doppler_cut_axis.step(doppler_khz, map_db[peak_row], where="mid", lw=1)
    doppler_cut_axis.axvline(TARGET_DOPPLER_HZ / 1e3, color="gray", lw=0.8, ls="--",
                             label="Target placed")
    doppler_cut_axis.set_xlabel("Doppler shift (kHz)")
    doppler_cut_axis.set_title(f"Doppler cut: row at {range_km[peak_row]:.2f} km")

    range_cut_axis.step(range_km, map_db[:, peak_column], where="mid", lw=1)
    range_cut_axis.axvline(range_km[target_delay], color="gray", lw=0.8, ls="--",
                           label="Target placed")
    range_cut_axis.set_xlabel("Range (km)")
    range_cut_axis.set_title(f"Range cut: column at {doppler_khz[peak_column]:.2f} kHz")

    for cut_axis in (doppler_cut_axis, range_cut_axis):
        cut_axis.set_ylim(-dynamic_range_db - 10, 3)
        cut_axis.set_ylabel("|CAF| (dB relative to peak)")
        cut_axis.grid(alpha=0.3)
        cut_axis.legend(loc="upper right", fontsize=8)
    cuts_figure.tight_layout()
    cuts_figure.savefig(cuts_path, dpi=150)

    print(
        f"Target placed at  {range_km[target_delay]:.2f} km, "
        f"{TARGET_DOPPLER_HZ / 1e3:.2f} kHz"
    )
    print(f"CAF target found  {range_km[peak_row]:.2f} km, {doppler_khz[peak_column]:.2f} kHz")
    print(f"Clutter placed at {', '.join(f'{range_km[d]:.2f}' for d in clutter_delays)} km, 0 Hz")
    print(f"Saved plots to {output_path} and {cuts_path}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CAF CPU baseline")
    parser.add_argument("--no-plot", action="store_true", help="skip the range-Doppler plot")
    parser.add_argument("--output", default="range_doppler_map.png", help="map plot file name; the cuts go to <name>_cuts.png")
    arguments = parser.parse_args()

    random_generator = np.random.default_rng(seed=0)
    if not arguments.no_plot:
        plot_range_doppler_map(random_generator, arguments.output)
    run_benchmark(random_generator)
