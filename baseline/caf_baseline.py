"""
CPU baseline for the Cross Ambiguity Function (CAF).

This is the software reference model for the FPGA accelerator. It follows the
structure of the Python prototype in Figure 4 of:

    C. D. Long, T. Brant, L. Brand, J. D. Ouellette, "Hardware-in-the-Loop
    Simulation with PYNQ for Radar Applications," Proc. 15th GNU Radio
    Conference, 2025.

The CAF from Eq. 1 of the paper is

    X(v) = IFFT( conj(FFT(x)) * FFT(y * exp(2j*pi*v*m / Fs)) ) / N

where x is the transmit signal, y is the receive signal, m is the sample
index, Fs is the sample rate, N is the number of samples, and v sweeps the
Doppler bins. The same three steps appear here as three functions:

    1. calculate_shifted_signals  -> frequency shifted copies of the receive signal
    2. cross_correlation          -> FFT based cross correlation
    3. calculate_cross_ambiguity  -> one cross correlation per Doppler bin
"""

import numpy as np


def calculate_shifted_signals(signal, num_of_samples, num_of_freq_bins, freq_step, sample_rate):
    """Build one frequency shifted copy of the signal for every Doppler bin.

    Row j holds the signal shifted by (j - num_of_freq_bins/2) * freq_step Hz,
    so the shifts are centered on 0 Hz.
    """
    shifted_signals = np.zeros((num_of_freq_bins, num_of_samples), dtype=np.complex64)
    sample_times = np.arange(num_of_samples) / sample_rate

    for freq_bin in range(num_of_freq_bins):
        freq_shift = (freq_bin - (num_of_freq_bins // 2)) * freq_step
        phase_shift = 2 * np.pi * freq_shift * sample_times
        shifted_signals[freq_bin] = signal * np.exp(1j * phase_shift)

    return shifted_signals


def cross_correlation(transmit_signal, receive_signal):
    """Circular cross correlation of the two signals, computed with FFTs.

    A copy of the transmit signal that shows up in the receive signal delayed
    by d samples produces a peak at index d.
    """
    transmit_fft = np.fft.fft(transmit_signal)
    receive_fft = np.fft.fft(receive_signal)
    correlation_fft = np.conj(transmit_fft) * receive_fft
    return np.fft.ifft(correlation_fft)


def calculate_cross_ambiguity(
    transmit_signal, receive_signal, num_of_samples, num_of_freq_bins, freq_step, sample_rate
):
    """Compute the range-Doppler map for a transmit / receive signal pair.

    Returns a real array of shape (num_of_samples, num_of_freq_bins). Rows are
    delay (range) bins and columns are Doppler bins. Use delay_axis() and
    doppler_axis() to label them.

    As in the paper, the output is the real part of the correlation. The
    complex result stays available from calculate_cross_ambiguity_complex()
    when the magnitude is needed instead.
    """
    cross_ambiguity = calculate_cross_ambiguity_complex(
        transmit_signal, receive_signal, num_of_samples, num_of_freq_bins, freq_step, sample_rate
    )
    return np.fliplr(cross_ambiguity.real.T)


def calculate_cross_ambiguity_complex(
    transmit_signal, receive_signal, num_of_samples, num_of_freq_bins, freq_step, sample_rate
):
    """Complex CAF of shape (num_of_freq_bins, num_of_samples), before the
    real / transpose / flip that produces the final range-Doppler map."""
    cross_ambiguity = np.zeros((num_of_freq_bins, num_of_samples), dtype=np.complex64)
    shifted_signals = calculate_shifted_signals(
        receive_signal, num_of_samples, num_of_freq_bins, freq_step, sample_rate
    )

    for freq_bin in range(num_of_freq_bins):
        cross_correlation_result = cross_correlation(transmit_signal, shifted_signals[freq_bin])
        cross_ambiguity[freq_bin] = cross_correlation_result / num_of_samples

    return cross_ambiguity


def doppler_axis(num_of_freq_bins, freq_step):
    """Target Doppler shift in Hz for each column of the range-Doppler map.

    Shifting the receive signal by -f cancels a target Doppler of +f, and the
    map columns are flipped, so column k corresponds to a target Doppler of
    (k - num_of_freq_bins/2 + 1) * freq_step.
    """
    return (np.arange(num_of_freq_bins) - (num_of_freq_bins // 2) + 1) * freq_step


def delay_axis(num_of_samples, sample_rate):
    """Time delay in seconds for each row of the range-Doppler map."""
    return np.arange(num_of_samples) / sample_rate
