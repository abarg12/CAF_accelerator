"""
Run the CAF CPU baseline: plot a range-Doppler map and time the computation.

    python baseline/run_baseline.py              # plot + benchmark
    python baseline/run_baseline.py --no-plot    # benchmark only

The benchmark sweeps square range-Doppler maps from 8x8 to 256x256, the same
sizes reported in the paper, so the numbers can be compared against the FPGA
core later.
"""

import argparse
import time

import numpy as np

from caf_baseline import calculate_cross_ambiguity, delay_axis, doppler_axis
from signal_generation import (
    delay_to_range_meters,
    generate_receive_signal,
    generate_transmit_signal,
)

SAMPLE_RATE = 1_000_000  # Hz
SIGNAL_TO_NOISE_RATIO_DB = 0
BENCHMARK_BIN_SIZES = [8, 16, 32, 64, 128, 256]
BENCHMARK_REPEATS = 20


def build_scenario(num_of_samples, num_of_freq_bins, random_generator):
    """Simulated transmit / receive pair with one target placed on the grid."""
    freq_step = SAMPLE_RATE / num_of_samples
    target_delay_in_samples = num_of_samples // 4
    target_doppler_hz = (num_of_freq_bins // 8) * freq_step

    transmit_signal = generate_transmit_signal(num_of_samples, random_generator)
    receive_signal = generate_receive_signal(
        transmit_signal,
        target_delay_in_samples,
        target_doppler_hz,
        SAMPLE_RATE,
        SIGNAL_TO_NOISE_RATIO_DB,
        random_generator,
    )
    return transmit_signal, receive_signal, freq_step, target_delay_in_samples, target_doppler_hz


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


def plot_range_doppler_map(random_generator, output_path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    num_of_samples = 256
    num_of_freq_bins = 256
    transmit_signal, receive_signal, freq_step, target_delay, target_doppler_hz = build_scenario(
        num_of_samples, num_of_freq_bins, random_generator
    )
    range_doppler_map = calculate_cross_ambiguity(
        transmit_signal, receive_signal, num_of_samples, num_of_freq_bins, freq_step, SAMPLE_RATE
    )

    doppler_hz = doppler_axis(num_of_freq_bins, freq_step)
    range_km = delay_to_range_meters(delay_axis(num_of_samples, SAMPLE_RATE)) / 1e3

    plt.figure(figsize=(7, 5))
    plt.pcolormesh(doppler_hz / 1e3, range_km, range_doppler_map, shading="auto")
    plt.xlabel("Doppler shift (kHz)")
    plt.ylabel("Range (km)")
    plt.title("CAF CPU baseline: range-Doppler map")
    plt.colorbar(label="Correlation (real part)")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)

    peak_row, peak_column = np.unravel_index(np.argmax(range_doppler_map), range_doppler_map.shape)
    print(f"Target placed at  {range_km[target_delay]:.2f} km, {target_doppler_hz / 1e3:.2f} kHz")
    print(f"CAF peak found at {range_km[peak_row]:.2f} km, {doppler_hz[peak_column] / 1e3:.2f} kHz")
    print(f"Saved plot to {output_path}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CAF CPU baseline")
    parser.add_argument("--no-plot", action="store_true", help="skip the range-Doppler plot")
    parser.add_argument("--output", default="range_doppler_map.png", help="plot file name")
    arguments = parser.parse_args()

    random_generator = np.random.default_rng(seed=0)
    if not arguments.no_plot:
        plot_range_doppler_map(random_generator, arguments.output)
    run_benchmark(random_generator)
