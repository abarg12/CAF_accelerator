/*
 * Cross Ambiguity Function (CAF), optimized Vitis HLS design.
 *
 * Computes the same thing as hls/baseline/caf.cpp,
 *
 *     X(v) = IFFT( conj(FFT(x)) * FFT(y * exp(2j*pi*v*m / Fs)) ) / N,
 *
 * with the same interface, but as a dataflow pipeline of five processes
 * connected by streams. Each process loops over the Doppler bins, so while
 * one bin is in the inverse FFT the next is in the forward FFT and the one
 * after that is being frequency shifted:
 *
 *   shift_signals -> fft_frames -> multiply_spectra -> fft_frames -> write_map
 *                    (forward)                         (inverse)
 *
 * The transmit FFT reuses the forward FFT: shift_signals sends the transmit
 * signal as frame 0, before the shifted receive signals, and
 * multiply_spectra keeps that first spectrum as its reference. This saves
 * a third FFT core.
 *
 * Every loop is pipelined at II=1. The FFT optimizations are described in
 * fft.cpp.
 */

#include <cmath>
#include "caf.h"

/*
 * Frame 0 is the transmit signal. Frames 1 .. NUM_OF_FREQ_BINS are the
 * receive signal multiplied by exp(2j*pi*freq_shift*m / sample_rate) for
 * each Doppler bin.
 */
static void shift_signals(
    DTYPE transmit_real[NUM_OF_SAMPLES],
    DTYPE transmit_imag[NUM_OF_SAMPLES],
    DTYPE receive_real[NUM_OF_SAMPLES],
    DTYPE receive_imag[NUM_OF_SAMPLES],
    DTYPE freq_step,
    DTYPE sample_rate,
    hls::stream<complex_sample> &output) {
transmit_loop:
	for (int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
#pragma HLS pipeline II=1
		complex_sample value;
		value.real = transmit_real[sample];
		value.imag = transmit_imag[sample];
		output.write(value);
	}

freq_bin_loop:
	for (int freq_bin = 0; freq_bin < NUM_OF_FREQ_BINS; freq_bin++) {
		DTYPE freq_shift = (freq_bin - NUM_OF_FREQ_BINS / 2) * freq_step;
		DTYPE phase_step = 6.283185307178f * freq_shift / sample_rate;

	shift_loop:
		for (int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
#pragma HLS pipeline II=1
			DTYPE phase_shift = phase_step * sample;
			DTYPE c = std::cos(phase_shift);
			DTYPE s = std::sin(phase_shift);
			complex_sample value;
			value.real = receive_real[sample] * c - receive_imag[sample] * s;
			value.imag = receive_real[sample] * s + receive_imag[sample] * c;
			output.write(value);
		}
	}
}

/*
 * Keeps the first spectrum (the transmit FFT) and multiplies each later
 * spectrum by its conjugate: conj(FFT(transmit)) * FFT(shifted receive).
 */
static void multiply_spectra(
    hls::stream<complex_sample> &input,
    hls::stream<complex_sample> &output) {
	DTYPE transmit_fft_real[NUM_OF_SAMPLES];
	DTYPE transmit_fft_imag[NUM_OF_SAMPLES];

reference_loop:
	for (int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
#pragma HLS pipeline II=1
		complex_sample value = input.read();
		transmit_fft_real[sample] = value.real;
		transmit_fft_imag[sample] = value.imag;
	}

freq_bin_loop:
	for (int freq_bin = 0; freq_bin < NUM_OF_FREQ_BINS; freq_bin++) {
	multiply_loop:
		for (int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
#pragma HLS pipeline II=1
			complex_sample value = input.read();
			complex_sample product;
			product.real = transmit_fft_real[sample] * value.real
			             + transmit_fft_imag[sample] * value.imag;
			product.imag = transmit_fft_real[sample] * value.imag
			             - transmit_fft_imag[sample] * value.real;
			output.write(product);
		}
	}
}

/*
 * Keeps the real part, scaled by 1/N for the inverse FFT and 1/N for the
 * CAF. The Doppler axis is flipped, as in np.fliplr().
 */
static void write_map(
    hls::stream<complex_sample> &input,
    DTYPE range_doppler_map[NUM_OF_SAMPLES * NUM_OF_FREQ_BINS]) {
	// 1/N^2 is a power of two, so this multiply is exact, like the
	// baseline's divide, but much cheaper
	const DTYPE scale = 1.0f / (NUM_OF_SAMPLES * NUM_OF_SAMPLES);

freq_bin_loop:
	for (int freq_bin = 0; freq_bin < NUM_OF_FREQ_BINS; freq_bin++) {
		int column = NUM_OF_FREQ_BINS - 1 - freq_bin;
	output_loop:
		for (int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
#pragma HLS pipeline II=1
			range_doppler_map[sample * NUM_OF_FREQ_BINS + column] = input.read().real * scale;
		}
	}
}

void cross_ambiguity(
    DTYPE transmit_real[NUM_OF_SAMPLES],
    DTYPE transmit_imag[NUM_OF_SAMPLES],
    DTYPE receive_real[NUM_OF_SAMPLES],
    DTYPE receive_imag[NUM_OF_SAMPLES],
    DTYPE freq_step,
    DTYPE sample_rate,
    DTYPE range_doppler_map[NUM_OF_SAMPLES * NUM_OF_FREQ_BINS]) {
#pragma HLS dataflow
	hls::stream<complex_sample> shifted("shifted");
	hls::stream<complex_sample> spectra("spectra");
	hls::stream<complex_sample> products("products");
	hls::stream<complex_sample> correlations("correlations");

	shift_signals(transmit_real, transmit_imag, receive_real, receive_imag,
	              freq_step, sample_rate, shifted);
	fft_frames(shifted, spectra, NUM_OF_FREQ_BINS + 1, false);
	multiply_spectra(spectra, products);
	fft_frames(products, correlations, NUM_OF_FREQ_BINS, true);
	write_map(correlations, range_doppler_map);
}
