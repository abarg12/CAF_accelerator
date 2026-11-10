/*
 * Self checking testbench for the CAF HLS baseline.
 *
 * It builds a noise-like transmit signal, makes a receive signal that is a
 * delayed and Doppler shifted copy with noise added, and runs the CAF. Two
 * checks are made:
 *
 *   1. The peak of the range-Doppler map sits at the delay and Doppler that
 *      were put into the receive signal.
 *   2. Every output matches a reference CAF computed directly from the
 *      correlation sum in double precision (no FFT).
 *
 * Returns 0 on success, which is what Vitis HLS C simulation expects.
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "caf.h"

const double PI = 3.14159265358979323846;
const double SAMPLE_RATE = 1000000.0;
const double FREQ_STEP = SAMPLE_RATE / NUM_OF_SAMPLES;
const int TARGET_DELAY_IN_SAMPLES = 40;
const int TARGET_DOPPLER_IN_BINS = 12;
const double NOISE_AMPLITUDE = 0.1;
const double TOLERANCE = 1e-3;

/* Uniform random value in [-1, 1) */
double random_value() {
	return 2.0 * rand() / (RAND_MAX + 1.0) - 1.0;
}

int main() {
	static DTYPE transmit_real[NUM_OF_SAMPLES], transmit_imag[NUM_OF_SAMPLES];
	static DTYPE receive_real[NUM_OF_SAMPLES], receive_imag[NUM_OF_SAMPLES];
	static DTYPE range_doppler_map[NUM_OF_SAMPLES * NUM_OF_FREQ_BINS];

	srand(1);
	for (int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
		transmit_real[sample] = random_value();
		transmit_imag[sample] = random_value();
	}

	// Receive signal: circular delay, Doppler shift, additive noise
	double target_doppler_hz = TARGET_DOPPLER_IN_BINS * FREQ_STEP;
	for (int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
		int delayed = (sample - TARGET_DELAY_IN_SAMPLES + NUM_OF_SAMPLES) % NUM_OF_SAMPLES;
		double phase = 2 * PI * target_doppler_hz * sample / SAMPLE_RATE;
		receive_real[sample] = transmit_real[delayed] * cos(phase) - transmit_imag[delayed] * sin(phase)
		                     + NOISE_AMPLITUDE * random_value();
		receive_imag[sample] = transmit_real[delayed] * sin(phase) + transmit_imag[delayed] * cos(phase)
		                     + NOISE_AMPLITUDE * random_value();
	}

	cross_ambiguity(transmit_real, transmit_imag, receive_real, receive_imag,
	                FREQ_STEP, SAMPLE_RATE, range_doppler_map);

	// Check 1: peak location
	int peak_row = 0, peak_column = 0;
	for (int row = 0; row < NUM_OF_SAMPLES; row++) {
		for (int column = 0; column < NUM_OF_FREQ_BINS; column++) {
			if (range_doppler_map[row * NUM_OF_FREQ_BINS + column] >
			    range_doppler_map[peak_row * NUM_OF_FREQ_BINS + peak_column]) {
				peak_row = row;
				peak_column = column;
			}
		}
	}
	// Column k corresponds to a target Doppler of (k - NUM_OF_FREQ_BINS/2 + 1) bins
	int peak_doppler_in_bins = peak_column - NUM_OF_FREQ_BINS / 2 + 1;
	printf("Target placed at  delay %d samples, Doppler %d bins\n",
	       TARGET_DELAY_IN_SAMPLES, TARGET_DOPPLER_IN_BINS);
	printf("CAF peak found at delay %d samples, Doppler %d bins\n", peak_row, peak_doppler_in_bins);
	bool peak_is_correct = (peak_row == TARGET_DELAY_IN_SAMPLES)
	                    && (peak_doppler_in_bins == TARGET_DOPPLER_IN_BINS);

	// Check 2: compare against the direct correlation sum
	double max_error = 0.0;
	static double shifted_real[NUM_OF_SAMPLES], shifted_imag[NUM_OF_SAMPLES];
	for (int freq_bin = 0; freq_bin < NUM_OF_FREQ_BINS; freq_bin++) {
		double freq_shift = (freq_bin - NUM_OF_FREQ_BINS / 2) * FREQ_STEP;
		for (int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
			double phase = 2 * PI * freq_shift * sample / SAMPLE_RATE;
			shifted_real[sample] = receive_real[sample] * cos(phase) - receive_imag[sample] * sin(phase);
			shifted_imag[sample] = receive_real[sample] * sin(phase) + receive_imag[sample] * cos(phase);
		}
		int column = NUM_OF_FREQ_BINS - 1 - freq_bin;
		for (int delay = 0; delay < NUM_OF_SAMPLES; delay++) {
			double sum_real = 0.0;
			for (int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
				int shifted_index = (sample + delay) % NUM_OF_SAMPLES;
				// real part of conj(transmit) * shifted receive
				sum_real += transmit_real[sample] * shifted_real[shifted_index]
				          + transmit_imag[sample] * shifted_imag[shifted_index];
			}
			double expected = sum_real / NUM_OF_SAMPLES;
			double error = fabs(expected - range_doppler_map[delay * NUM_OF_FREQ_BINS + column]);
			if (error > max_error) {
				max_error = error;
			}
		}
	}
	printf("Largest difference from the reference CAF: %g (tolerance %g)\n", max_error, TOLERANCE);
	bool values_are_correct = max_error < TOLERANCE;

	if (peak_is_correct && values_are_correct) {
		printf("PASS\n");
		return 0;
	}
	printf("FAIL\n");
	return 1;
}
