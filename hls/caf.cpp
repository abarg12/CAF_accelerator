/*
 * Cross Ambiguity Function (CAF), unoptimized Vitis HLS baseline.
 *
 * This is a direct C++ port of the Python reference model in
 * baseline/caf_baseline.py, which follows Figure 4 of Long et al.,
 * "Hardware-in-the-Loop Simulation with PYNQ for Radar Applications"
 * (GNU Radio Conference 2025):
 *
 *     X(v) = IFFT( conj(FFT(x)) * FFT(y * exp(2j*pi*v*m / Fs)) ) / N
 *
 * There are no HLS pragmas here on purpose. This version gives us the
 * starting latency and resource numbers that later optimizations
 * (pipelining, dataflow, the Xilinx FFT IP) are measured against.
 *
 * The FFT in fft.cpp is adapted from https://github.com/KastnerRG/pp4fpgas
 */

#include <cmath>
#include "caf.h"

/* Multiply the signal by exp(2j*pi*freq_shift*m / sample_rate) */
void shift_signal(
    DTYPE signal_real[NUM_OF_SAMPLES],
    DTYPE signal_imag[NUM_OF_SAMPLES],
    DTYPE freq_shift,
    DTYPE sample_rate,
    DTYPE shifted_real[NUM_OF_SAMPLES],
    DTYPE shifted_imag[NUM_OF_SAMPLES]) {
	DTYPE phase_step = 6.283185307178 * freq_shift / sample_rate;

shift_loop:
	for (int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
		DTYPE phase_shift = phase_step * sample;
		DTYPE c = cos(phase_shift);
		DTYPE s = sin(phase_shift);
		shifted_real[sample] = signal_real[sample] * c - signal_imag[sample] * s;
		shifted_imag[sample] = signal_real[sample] * s + signal_imag[sample] * c;
	}
}

/*
 * Cross correlation of the transmit signal (already in the frequency domain)
 * with a receive signal: IFFT( conj(FFT(transmit)) * FFT(receive) ).
 * The result overwrites the receive arrays.
 */
void cross_correlation(
    DTYPE transmit_fft_real[NUM_OF_SAMPLES],
    DTYPE transmit_fft_imag[NUM_OF_SAMPLES],
    DTYPE receive_real[NUM_OF_SAMPLES],
    DTYPE receive_imag[NUM_OF_SAMPLES]) {
	fft(receive_real, receive_imag, false);

multiply_loop:
	for (int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
		// conj(transmit_fft) * receive_fft
		DTYPE product_real = transmit_fft_real[sample] * receive_real[sample]
		                   + transmit_fft_imag[sample] * receive_imag[sample];
		DTYPE product_imag = transmit_fft_real[sample] * receive_imag[sample]
		                   - transmit_fft_imag[sample] * receive_real[sample];
		receive_real[sample] = product_real;
		receive_imag[sample] = product_imag;
	}

	fft(receive_real, receive_imag, true);
}

void cross_ambiguity(
    DTYPE transmit_real[NUM_OF_SAMPLES],
    DTYPE transmit_imag[NUM_OF_SAMPLES],
    DTYPE receive_real[NUM_OF_SAMPLES],
    DTYPE receive_imag[NUM_OF_SAMPLES],
    DTYPE freq_step,
    DTYPE sample_rate,
    DTYPE range_doppler_map[NUM_OF_SAMPLES * NUM_OF_FREQ_BINS]) {
	DTYPE transmit_fft_real[NUM_OF_SAMPLES];
	DTYPE transmit_fft_imag[NUM_OF_SAMPLES];
	DTYPE shifted_real[NUM_OF_SAMPLES];
	DTYPE shifted_imag[NUM_OF_SAMPLES];

	// The transmit FFT is the same for every Doppler bin, so compute it once
copy_loop:
	for (int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
		transmit_fft_real[sample] = transmit_real[sample];
		transmit_fft_imag[sample] = transmit_imag[sample];
	}
	fft(transmit_fft_real, transmit_fft_imag, false);

freq_bin_loop:
	for (int freq_bin = 0; freq_bin < NUM_OF_FREQ_BINS; freq_bin++) {
		DTYPE freq_shift = (freq_bin - NUM_OF_FREQ_BINS / 2) * freq_step;

		shift_signal(receive_real, receive_imag, freq_shift, sample_rate, shifted_real, shifted_imag);
		cross_correlation(transmit_fft_real, transmit_fft_imag, shifted_real, shifted_imag);

		// Keep the real part, scaled by 1/N for the inverse FFT and 1/N for
		// the CAF. The Doppler axis is flipped, as in np.fliplr().
		int column = NUM_OF_FREQ_BINS - 1 - freq_bin;
	output_loop:
		for (int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
			range_doppler_map[sample * NUM_OF_FREQ_BINS + column] =
			    shifted_real[sample] / (NUM_OF_SAMPLES * NUM_OF_SAMPLES);
		}
	}
}
