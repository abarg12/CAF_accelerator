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

from caf_baseline import calculate_cross_ambiguity, calculate_cross_ambiguity_complex
from caf_baseline import delay_axis, doppler_axis
from signal_generation import (
    add_stationary_clutter,
    delay_to_range_meters,
    generate_fm_transmit_signal,
    generate_receive_signal,
    range_meters_to_delay_samples,
)

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


def build_scenario(num_of_samples, num_of_freq_bins, random_generator):
    """Simulated transmit / receive pair with one moving target and clutter.

    Delays wrap around for maps too small to hold the full range.
    """
    freq_step = DOPPLER_EXTENT_HZ / num_of_freq_bins
    target_delay = range_meters_to_delay_samples(TARGET_RANGE_METERS, SAMPLE_RATE) % num_of_samples
    clutter_delays = [
        range_meters_to_delay_samples(r, SAMPLE_RATE) % num_of_samples for r, _ in CLUTTER
    ]
    clutter_amplitudes = [amplitude for _, amplitude in CLUTTER]

    transmit_signal = generate_fm_transmit_signal(num_of_samples, SAMPLE_RATE, random_generator)
    receive_signal = generate_receive_signal(
        transmit_signal,
        target_delay,
        TARGET_DOPPLER_HZ,
        SAMPLE_RATE,
        SIGNAL_TO_NOISE_RATIO_DB,
        random_generator,
    )
    receive_signal = add_stationary_clutter(
        receive_signal, transmit_signal, clutter_delays, clutter_amplitudes
    )
    return transmit_signal, receive_signal, freq_step, target_delay, clutter_delays


def run_benchmark(random_generator):
    print(f"{'bins':>10} {'mean (ms)':>12} {'min (ms)':>12}")
    for bin_size in BENCHMARK_BIN_SIZES:
        transmit_signal, receive_signal, freq_step, _, _ = build_scenario(
            bin_size, bin_size, random_generator
        )
        run_times_ms = []
        for _ in range(BENCHMARK_REPEATS):
            start_time = time.perf_counter()
            calculate_cross_ambiguity(
                transmit_signal, receive_signal, bin_size, bin_size, freq_step, SAMPLE_RATE
            )
            run_times_ms.append((time.perf_counter() - start_time) * 1e3)
        print(
            f"{f'{bin_size}x{bin_size}':>10} "
            f"{np.mean(run_times_ms):>12.3f} {np.min(run_times_ms):>12.3f}"
        )


def range_doppler_magnitude_db(transmit_signal, receive_signal, num_of_samples, num_of_freq_bins,
                               freq_step):
    """|CAF| in dB relative to the map's peak, in the same orientation as
    calculate_cross_ambiguity(): rows are range bins, columns are Doppler bins.

    The paper and the FPGA core output the real part, which swings with the
    phase of each return. The magnitude does not, so it is easier to read.
    """
    cross_ambiguity = calculate_cross_ambiguity_complex(
        transmit_signal, receive_signal, num_of_samples, num_of_freq_bins, freq_step, SAMPLE_RATE
    )
    magnitude = np.fliplr(np.abs(cross_ambiguity).T)
    return 20 * np.log10(magnitude / magnitude.max() + 1e-12)


def cuts_output_path(map_output_path):
    """File name for the cuts figure: range_doppler_map.png -> range_doppler_map_cuts.png"""
    stem, extension = os.path.splitext(map_output_path)
    return f"{stem}_cuts{extension or '.png'}"


def plot_range_doppler_map(random_generator, output_path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    num_of_samples = 256
    num_of_freq_bins = 256
    dynamic_range_db = 20

    transmit_signal, receive_signal, freq_step, target_delay, clutter_delays = build_scenario(
        num_of_samples, num_of_freq_bins, random_generator
    )
    map_db = range_doppler_magnitude_db(
        transmit_signal, receive_signal, num_of_samples, num_of_freq_bins, freq_step
    )
    doppler_khz = doppler_axis(num_of_freq_bins, freq_step) / 1e3
    range_km = delay_to_range_meters(delay_axis(num_of_samples, SAMPLE_RATE)) / 1e3

    # Detect the moving target as the strongest return outside the zero
    # Doppler clutter band (one Doppler resolution cell on each side of 0 Hz).
    doppler_resolution_khz = SAMPLE_RATE / num_of_samples / 1e3
    clutter_band = np.abs(doppler_khz) < doppler_resolution_khz
    moving_map_db = np.where(clutter_band[np.newaxis, :], -np.inf, map_db)
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
    cuts_path = cuts_output_path(output_path)
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
