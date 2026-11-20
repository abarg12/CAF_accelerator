"""Checks that the CAF baseline recovers known targets and clutter."""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "baseline"))

import run_baseline  # noqa: E402
from caf_baseline import cross_ambiguity  # noqa: E402

SAMPLE_RATE = 1_000_000


def white_noise(num_of_samples, random_generator):
    """Unit power complex white noise. Its ambiguity function is a one bin
    "thumbtack", which makes exact delay and Doppler checks easy."""
    return (
        (random_generator.standard_normal(num_of_samples)
         + 1j * random_generator.standard_normal(num_of_samples)) / np.sqrt(2)
    ).astype(np.complex64)


@pytest.mark.parametrize(
    "num_of_samples, delay, doppler_bins", [(256, 40, 12), (256, 100, -30), (128, 5, 0)]
)
def test_recovers_target_delay_and_doppler(num_of_samples, delay, doppler_bins):
    random_generator = np.random.default_rng(seed=1)
    freq_step = SAMPLE_RATE / num_of_samples
    sample_times = np.arange(num_of_samples) / SAMPLE_RATE

    transmit_signal = white_noise(num_of_samples, random_generator)
    receive_signal = np.roll(transmit_signal, delay) * np.exp(
        2j * np.pi * doppler_bins * freq_step * sample_times
    ) + 0.3 * white_noise(num_of_samples, random_generator)

    caf = cross_ambiguity(transmit_signal, receive_signal, num_of_samples, freq_step, SAMPLE_RATE)
    peak_row, peak_column = np.unravel_index(np.argmax(caf.real), caf.shape)
    assert peak_row == delay
    assert peak_column - num_of_samples // 2 + 1 == doppler_bins


def test_output_shape_is_samples_by_freq_bins():
    signal = white_noise(256, np.random.default_rng(seed=1))
    assert cross_ambiguity(signal, signal, 64, 1000, SAMPLE_RATE).shape == (256, 64)


def test_zero_doppler_column_matches_direct_correlation_sum():
    random_generator = np.random.default_rng(seed=2)
    num_of_samples = 32
    transmit_signal = white_noise(num_of_samples, random_generator)
    receive_signal = white_noise(num_of_samples, random_generator)

    direct = np.array(
        [
            np.sum(np.conj(transmit_signal) * np.roll(receive_signal, -lag))
            for lag in range(num_of_samples)
        ]
    ) / num_of_samples
    caf = cross_ambiguity(transmit_signal, receive_signal, 4, 1000, SAMPLE_RATE)
    # Column k is a Doppler of (k - 4/2 + 1) bins, so column 1 is 0 Hz
    assert np.allclose(caf[:, 1], direct, atol=1e-4)


def test_stationary_clutter_appears_at_zero_doppler():
    num_of_samples = 256
    transmit_signal = white_noise(num_of_samples, np.random.default_rng(seed=4))
    receive_signal = np.roll(transmit_signal, 30) + 0.5 * np.roll(transmit_signal, 90)

    caf = cross_ambiguity(
        transmit_signal, receive_signal, num_of_samples, SAMPLE_RATE / num_of_samples, SAMPLE_RATE
    )
    zero_doppler_column = num_of_samples // 2 - 1
    assert np.isclose(caf.real[30, zero_doppler_column], 1.0, atol=0.05)
    assert np.isclose(caf.real[90, zero_doppler_column], 0.5, atol=0.05)


def test_fm_transmit_signal_is_constant_envelope_and_band_limited():
    num_of_samples = 4096
    signal, _, _, _ = run_baseline.make_signals(num_of_samples, np.random.default_rng(seed=3))
    assert np.allclose(np.abs(signal), 1, atol=1e-5)

    # Carson's rule: about 2 * (75 + 53) = 256 kHz. Almost all power lies within
    # +/-200 kHz, far less than the +/-500 kHz a white noise signal would cover.
    power = np.abs(np.fft.fft(signal)) ** 2
    frequencies = np.fft.fftfreq(num_of_samples, d=1 / SAMPLE_RATE)
    assert power[np.abs(frequencies) <= 200_000].sum() / power.sum() > 0.99


def test_paper_scenario_recovers_moving_target():
    num_of_samples = 256
    freq_step = run_baseline.DOPPLER_EXTENT_HZ / num_of_samples
    transmit_signal, receive_signal, target_delay, _ = run_baseline.make_signals(
        num_of_samples, np.random.default_rng(seed=0)
    )
    magnitude = np.abs(
        cross_ambiguity(transmit_signal, receive_signal, num_of_samples, freq_step, SAMPLE_RATE)
    )
    doppler_hz = (np.arange(num_of_samples) - num_of_samples // 2 + 1) * freq_step
    magnitude[:, np.abs(doppler_hz) < SAMPLE_RATE / num_of_samples] = 0  # ignore clutter
    peak_row, peak_column = np.unravel_index(np.argmax(magnitude), magnitude.shape)

    assert abs(peak_row - target_delay) <= 1
    # Doppler resolution is sample_rate / num_of_samples = 3.9 kHz
    assert abs(doppler_hz[peak_column] - run_baseline.TARGET_DOPPLER_HZ) < 1_000
