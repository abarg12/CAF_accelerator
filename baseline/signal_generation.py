"""
Simulated transmit and receive signals for testing the CAF baseline.

Two transmit signals are available:

    generate_transmit_signal     complex white noise. Its ambiguity function is
                                 a one-bin "thumbtack", which makes it a strict
                                 test signal for checking the CAF math.
    generate_fm_transmit_signal  a simulated FM broadcast signal. The paper
                                 records an FM radio station with an SDR and
                                 uses it as the transmit signal, so this is the
                                 one to use for paper-like range-Doppler maps.

The receive signal is a copy of the transmit signal with the three target
effects described in the paper: a time delay (range), a Doppler frequency shift
(velocity), and additive noise. add_stationary_clutter() adds the zero-Doppler
returns from stationary objects that the paper simulates in Section 2.6.
"""

import numpy as np

SPEED_OF_LIGHT = 299_792_458.0  # m/s


def generate_transmit_signal(num_of_samples, random_generator):
    """Unit power complex white noise."""
    real_part = random_generator.standard_normal(num_of_samples)
    imag_part = random_generator.standard_normal(num_of_samples)
    return ((real_part + 1j * imag_part) / np.sqrt(2)).astype(np.complex64)


def generate_fm_transmit_signal(
    num_of_samples,
    sample_rate,
    random_generator,
    message_bandwidth_hz=53_000,
    frequency_deviation_hz=75_000,
):
    """Unit power FM broadcast signal at complex baseband.

    The message stands in for the FM stereo multiplex (MPX) baseband, which
    carries the audio, stereo pilot, and stereo subcarrier out to about 53 kHz.
    It is white noise low pass filtered to message_bandwidth_hz,
    scaled to a peak of 1, and frequency modulated with the broadcast FM
    deviation of 75 kHz. By Carson's rule the signal occupies about
    2 * (75 + 53) = 256 kHz, so the range peak spreads over a few bins instead
    of the single bin a white noise transmit signal gives.
    """
    white_noise = random_generator.standard_normal(num_of_samples)
    noise_spectrum = np.fft.fft(white_noise)
    frequencies = np.fft.fftfreq(num_of_samples, d=1 / sample_rate)
    noise_spectrum[np.abs(frequencies) > message_bandwidth_hz] = 0
    message = np.fft.ifft(noise_spectrum).real
    message /= np.max(np.abs(message))

    phase = 2 * np.pi * frequency_deviation_hz * np.cumsum(message) / sample_rate
    return np.exp(1j * phase).astype(np.complex64)


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


def add_stationary_clutter(receive_signal, transmit_signal, delays_in_samples, amplitudes):
    """Add zero Doppler returns from stationary objects at the given delays.

    Each clutter return is a circularly delayed, scaled copy of the transmit
    signal with no frequency shift, so it shows up in the range-Doppler map at
    0 Hz.
    """
    clutter = np.zeros(len(transmit_signal), dtype=np.complex128)
    for delay_in_samples, amplitude in zip(delays_in_samples, amplitudes):
        clutter += amplitude * np.roll(transmit_signal, delay_in_samples)
    return (receive_signal + clutter).astype(np.complex64)


def range_meters_to_delay_samples(range_meters, sample_rate):
    """Nearest whole sample delay for a monostatic range."""
    return int(round(2 * range_meters / SPEED_OF_LIGHT * sample_rate))


def delay_to_range_meters(delay_in_seconds):
    """Monostatic range for a round trip delay."""
    return SPEED_OF_LIGHT * delay_in_seconds / 2
