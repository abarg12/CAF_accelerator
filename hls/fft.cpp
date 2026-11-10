/*
 * Software style radix-2 FFT, used as the unoptimized HLS baseline.
 *
 * Adapted from the fft_sw example (examples/fft_sw.cpp) in
 * "Parallel Programming for FPGAs" by R. Kastner, J. Matai, and S. Neuendorffer:
 *
 *     https://github.com/KastnerRG/pp4fpgas
 *
 * The original is licensed under Creative Commons Attribution 4.0
 * International (CC BY 4.0), https://creativecommons.org/licenses/by/4.0/
 *
 * Changes from the original: the size constants were renamed, and an
 * `inverse` argument was added so the same function also computes the
 * inverse FFT.
 */

#include <cmath>
#include "caf.h"

unsigned int reverse_bits(unsigned int input) {
	int i, rev = 0;
	for (i = 0; i < LOG2_NUM_OF_SAMPLES; i++) {
		rev = (rev << 1) | (input & 1);
		input = input >> 1;
	}
	return rev;
}

void bit_reverse(DTYPE X_R[NUM_OF_SAMPLES], DTYPE X_I[NUM_OF_SAMPLES]) {
	unsigned int reversed;
	unsigned int i;
	DTYPE temp;

	for (i = 0; i < NUM_OF_SAMPLES; i++) {
		reversed = reverse_bits(i); // Find the bit reversed index
		if (i <= reversed) {
			// Swap the real values
			temp = X_R[i];
			X_R[i] = X_R[reversed];
			X_R[reversed] = temp;

			// Swap the imaginary values
			temp = X_I[i];
			X_I[i] = X_I[reversed];
			X_I[reversed] = temp;
		}
	}
}

void fft(DTYPE X_R[NUM_OF_SAMPLES], DTYPE X_I[NUM_OF_SAMPLES], bool inverse) {
	DTYPE temp_R; // temporary storage complex variable
	DTYPE temp_I; // temporary storage complex variable
	int i, j;     // loop indexes
	int i_lower;  // Index of lower point in butterfly
	int stage, DFTpts;
	int numBF;    // Butterfly Width
	DTYPE a, e, c, s;

	// The inverse FFT uses the opposite sign in the twiddle factor exponent
	DTYPE two_pi = inverse ? 6.283185307178 : -6.283185307178;

	bit_reverse(X_R, X_I);

stage_loop:
	for (stage = 1; stage <= LOG2_NUM_OF_SAMPLES; stage++) { // Do log2(N) stages of butterflies
		DFTpts = 1 << stage;  // DFT = 2^stage = points in sub DFT
		numBF = DFTpts / 2;   // Butterfly WIDTHS in sub-DFT
		e = two_pi / DFTpts;
		a = 0.0;
	// Perform butterflies for j-th stage
	butterfly_loop:
		for (j = 0; j < numBF; j++) {
			c = cos(a);
			s = sin(a);
			a = a + e;
		// Compute butterflies that use same W**k
		dft_loop:
			for (i = j; i < NUM_OF_SAMPLES; i += DFTpts) {
				i_lower = i + numBF; // index of lower point in butterfly
				temp_R = X_R[i_lower] * c - X_I[i_lower] * s;
				temp_I = X_I[i_lower] * c + X_R[i_lower] * s;
				X_R[i_lower] = X_R[i] - temp_R;
				X_I[i_lower] = X_I[i] - temp_I;
				X_R[i] = X_R[i] + temp_R;
				X_I[i] = X_I[i] + temp_I;
			}
		}
	}
}
