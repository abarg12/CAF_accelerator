"""
CPU baseline for the Cross Ambiguity Function (CAF).

This is the software reference model for the FPGA accelerator. It follows the
Python prototype in Figure 4 of:

    C. D. Long, T. Brant, L. Brand, J. D. Ouellette, "Hardware-in-the-Loop
    Simulation with PYNQ for Radar Applications," Proc. 15th GNU Radio
    Conference, 2025.

The CAF from Eq. 1 of the paper is

    X(v) = IFFT( conj(FFT(x)) * FFT(y * exp(2j*pi*v*m / Fs)) ) / N

where x is the transmit signal, y is the receive signal, m is the sample
index, Fs is the sample rate, N is the number of samples, and v sweeps the
Doppler bins.
"""

import numpy as np


def cross_ambiguity(transmit_signal, receive_signal, num_of_freq_bins, freq_step, sample_rate):
    """Complex range-Doppler map of shape (num_of_samples, num_of_freq_bins).

    Rows are delay bins: row d is a delay of d / sample_rate seconds.
    Columns are Doppler bins: column k is a target Doppler of
    (k - num_of_freq_bins/2 + 1) * freq_step Hz.

    The FPGA core outputs the real part, as in the paper.
    """
    num_of_samples = len(transmit_signal)
    sample_times = np.arange(num_of_samples) / sample_rate
    transmit_fft_conj = np.conj(np.fft.fft(transmit_signal))

    caf = np.zeros((num_of_samples, num_of_freq_bins), dtype=np.complex64)
    for freq_bin in range(num_of_freq_bins):
        # Shifting the receive signal by -f cancels a target Doppler of +f.
        # Columns are filled in reverse, as np.fliplr() does in the paper.
        freq_shift = (freq_bin - num_of_freq_bins // 2) * freq_step
        shifted = receive_signal * np.exp(2j * np.pi * freq_shift * sample_times)
        correlation = np.fft.ifft(transmit_fft_conj * np.fft.fft(shifted))
        caf[:, num_of_freq_bins - 1 - freq_bin] = correlation / num_of_samples
    return caf
