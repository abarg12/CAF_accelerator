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
