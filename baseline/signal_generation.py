"""
Simulated transmit and receive signals for testing the CAF baseline.

The transmit signal is complex white noise, which stands in for the noise-like
illuminator used by a continuous wave radar. The receive signal is a copy of
the transmit signal with the three target effects described in the paper:
a time delay (range), a Doppler frequency shift (velocity), and additive noise.
"""

import numpy as np

SPEED_OF_LIGHT = 299_792_458.0  # m/s


def generate_transmit_signal(num_of_samples, random_generator):
    """Unit power complex white noise."""
    real_part = random_generator.standard_normal(num_of_samples)
    imag_part = random_generator.standard_normal(num_of_samples)
    return ((real_part + 1j * imag_part) / np.sqrt(2)).astype(np.complex64)


def generate_receive_signal(
    transmit_signal,
    delay_in_samples,
    doppler_shift_hz,
    sample_rate,
    signal_to_noise_ratio_db,
    random_generator,
):
    """Delay, Doppler shift, and add noise to the transmit signal.

    The delay is circular so that it lines up with the circular correlation
    that the FFT based CAF computes.
    """
    num_of_samples = len(transmit_signal)
    sample_times = np.arange(num_of_samples) / sample_rate

    delayed_signal = np.roll(transmit_signal, delay_in_samples)
    target_return = delayed_signal * np.exp(2j * np.pi * doppler_shift_hz * sample_times)

    noise_power = 10 ** (-signal_to_noise_ratio_db / 10)
    noise = np.sqrt(noise_power / 2) * (
        random_generator.standard_normal(num_of_samples)
        + 1j * random_generator.standard_normal(num_of_samples)
    )
    return (target_return + noise).astype(np.complex64)


def delay_to_range_meters(delay_in_seconds):
    """Monostatic range for a round trip delay."""
    return SPEED_OF_LIGHT * delay_in_seconds / 2
