"""Checks that the CAF baseline recovers a known target delay and Doppler."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "baseline"))

from caf_baseline import (  # noqa: E402
    calculate_cross_ambiguity,
    cross_correlation,
    delay_axis,
    doppler_axis,
)
from signal_generation import generate_receive_signal, generate_transmit_signal  # noqa: E402

SAMPLE_RATE = 1_000_000


def find_peak(num_of_samples, num_of_freq_bins, delay_in_samples, doppler_bins):
    random_generator = np.random.default_rng(seed=1)
    freq_step = SAMPLE_RATE / num_of_samples
    doppler_shift_hz = doppler_bins * freq_step

    transmit_signal = generate_transmit_signal(num_of_samples, random_generator)
    receive_signal = generate_receive_signal(
        transmit_signal, delay_in_samples, doppler_shift_hz, SAMPLE_RATE, 10, random_generator
    )
    range_doppler_map = calculate_cross_ambiguity(
        transmit_signal, receive_signal, num_of_samples, num_of_freq_bins, freq_step, SAMPLE_RATE
    )
    peak_row, peak_column = np.unravel_index(np.argmax(range_doppler_map), range_doppler_map.shape)
    return (
        range_doppler_map.shape,
        delay_axis(num_of_samples, SAMPLE_RATE)[peak_row] * SAMPLE_RATE,
        doppler_axis(num_of_freq_bins, freq_step)[peak_column] / freq_step,
    )


def test_output_shape_is_samples_by_freq_bins():
    shape, _, _ = find_peak(256, 64, 10, 0)
    assert shape == (256, 64)


def test_recovers_delay_and_positive_doppler():
    _, delay, doppler = find_peak(256, 256, 40, 12)
    assert round(delay) == 40
    assert round(doppler) == 12


def test_recovers_delay_and_negative_doppler():
    _, delay, doppler = find_peak(256, 256, 100, -30)
    assert round(delay) == 100
    assert round(doppler) == -30


def test_recovers_stationary_target():
    _, delay, doppler = find_peak(128, 128, 5, 0)
    assert round(delay) == 5
    assert round(doppler) == 0


def test_cross_correlation_matches_direct_sum():
    random_generator = np.random.default_rng(seed=2)
    num_of_samples = 32
    transmit_signal = generate_transmit_signal(num_of_samples, random_generator)
    receive_signal = generate_transmit_signal(num_of_samples, random_generator)

    direct = np.array(
        [
            np.sum(np.conj(transmit_signal) * np.roll(receive_signal, -lag))
            for lag in range(num_of_samples)
        ]
    )
    assert np.allclose(cross_correlation(transmit_signal, receive_signal), direct, atol=1e-3)


def test_fm_transmit_signal_is_constant_envelope_and_band_limited():
    from signal_generation import generate_fm_transmit_signal

    num_of_samples = 4096
    signal = generate_fm_transmit_signal(
        num_of_samples, SAMPLE_RATE, np.random.default_rng(seed=3)
    )
    assert np.allclose(np.abs(signal), 1, atol=1e-5)

    # Carson's rule: about 2 * (75 + 53) = 256 kHz. Almost all power lies within
    # +/-200 kHz, far less than the +/-500 kHz a white noise signal would cover.
    power = np.abs(np.fft.fft(signal)) ** 2
    frequencies = np.fft.fftfreq(num_of_samples, d=1 / SAMPLE_RATE)
    assert power[np.abs(frequencies) <= 200_000].sum() / power.sum() > 0.99


def test_stationary_clutter_appears_at_zero_doppler():
    from signal_generation import add_stationary_clutter

    num_of_samples = 256
    freq_step = SAMPLE_RATE / num_of_samples
    transmit_signal = generate_transmit_signal(num_of_samples, np.random.default_rng(seed=4))
    receive_signal = add_stationary_clutter(
        np.zeros(num_of_samples, dtype=np.complex64), transmit_signal, [30, 90], [1.0, 0.5]
    )
    range_doppler_map = calculate_cross_ambiguity(
        transmit_signal, receive_signal, num_of_samples, num_of_samples, freq_step, SAMPLE_RATE
    )
    zero_doppler_column = int(np.argmin(np.abs(doppler_axis(num_of_samples, freq_step))))
    assert np.isclose(range_doppler_map[30, zero_doppler_column], 1.0, atol=0.05)
    assert np.isclose(range_doppler_map[90, zero_doppler_column], 0.5, atol=0.05)


def test_paper_scenario_recovers_moving_target():
    import run_baseline

    num_of_samples = 256
    transmit_signal, receive_signal, freq_step, target_delay, _ = run_baseline.build_scenario(
        num_of_samples, num_of_samples, np.random.default_rng(seed=0)
    )
    map_db = run_baseline.range_doppler_magnitude_db(
        transmit_signal, receive_signal, num_of_samples, num_of_samples, freq_step
    )
    doppler_hz = doppler_axis(num_of_samples, freq_step)
    map_db[:, np.abs(doppler_hz) < SAMPLE_RATE / num_of_samples] = -np.inf  # ignore clutter
    peak_row, peak_column = np.unravel_index(np.argmax(map_db), map_db.shape)

    assert abs(peak_row - target_delay) <= 1
    # Doppler resolution is sample_rate / num_of_samples = 3.9 kHz
    assert abs(doppler_hz[peak_column] - run_baseline.TARGET_DOPPLER_HZ) < 1_000
