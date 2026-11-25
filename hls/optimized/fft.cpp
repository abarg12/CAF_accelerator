/*
 * Streaming radix-2 FFT for the optimized CAF.
 *
 * Compared with the baseline FFT in hls/baseline/fft.cpp:
 *
 *   - Twiddle factors come from a ROM (twiddles.h) instead of std::cos and
 *     std::sin. In the baseline every butterfly group called both, about
 *     56 cycles each time.
 *   - Each stage is one loop over all N/2 butterflies, pipelined at II=1.
 *     The baseline's inner loop ran at II=25 because the tool could not
 *     tell that butterflies in one stage touch different samples; the
 *     dependence pragma below tells it so.
 *   - The samples are split into two memory banks by the parity of their
 *     index (XOR of its bits). The two samples of a butterfly differ in
 *     exactly one bit, so they are always in different banks, and each
 *     butterfly does one read and one write per bank. That fits the two
 *     ports of a block RAM, which is what allows II=1.
 *   - The bit reversal is done while the frame is loaded, by writing
 *     sample m to address reverse(m). Bit reversal keeps the number of one
 *     bits, so a sample stays in the same bank.
 *   - The inverse FFT uses IFFT(X) = conj(FFT(conj(X))), so both directions
 *     share the same twiddle ROM and butterfly.
 */

#include "caf.h"
#include "twiddles.h"

#define HALF_OF_SAMPLES (NUM_OF_SAMPLES / 2)

static unsigned int reverse_bits(unsigned int input) {
	unsigned int reversed = 0;
	for (int bit = 0; bit < LOG2_NUM_OF_SAMPLES; bit++) {
#pragma HLS unroll
		reversed = (reversed << 1) | (input & 1);
		input >>= 1;
	}
	return reversed;
}

/* Bank of a sample: XOR of all bits of its index */
static unsigned int bank_of(unsigned int index) {
	unsigned int parity = 0;
	for (int bit = 0; bit < LOG2_NUM_OF_SAMPLES; bit++) {
#pragma HLS unroll
		parity ^= (index >> bit) & 1;
	}
	return parity;
}

void fft_frames(
    hls::stream<complex_sample> &input,
    hls::stream<complex_sample> &output,
    int num_of_frames,
    bool inverse) {
	// Sample i lives in bank bank_of(i) at address i / 2
	DTYPE bank0_real[HALF_OF_SAMPLES], bank0_imag[HALF_OF_SAMPLES];
	DTYPE bank1_real[HALF_OF_SAMPLES], bank1_imag[HALF_OF_SAMPLES];

frame_loop:
	for (int frame = 0; frame < num_of_frames; frame++) {
#pragma HLS loop_tripcount min=NUM_OF_FREQ_BINS max=NUM_OF_FREQ_BINS+1

	load_loop:
		for (unsigned int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
#pragma HLS pipeline II=1
			complex_sample value = input.read();
			DTYPE imag = inverse ? -value.imag : value.imag;
			unsigned int index = reverse_bits(sample);
			if (bank_of(index)) {
				bank1_real[index / 2] = value.real;
				bank1_imag[index / 2] = imag;
			} else {
				bank0_real[index / 2] = value.real;
				bank0_imag[index / 2] = imag;
			}
		}

	stage_loop:
		for (int stage = 0; stage < LOG2_NUM_OF_SAMPLES; stage++) {
			// The pipeline must drain between stages, because a stage reads
			// what the previous one wrote. Flattening would remove that.
#pragma HLS loop_flatten off
			unsigned int span = 1u << stage; // distance between butterfly inputs

		butterfly_loop:
			for (unsigned int butterfly = 0; butterfly < HALF_OF_SAMPLES; butterfly++) {
#pragma HLS pipeline II=1
#pragma HLS dependence variable=bank0_real type=inter dependent=false
#pragma HLS dependence variable=bank0_imag type=inter dependent=false
#pragma HLS dependence variable=bank1_real type=inter dependent=false
#pragma HLS dependence variable=bank1_imag type=inter dependent=false
				unsigned int position = butterfly & (span - 1);           // position in its sub-DFT
				unsigned int upper = ((butterfly >> stage) << (stage + 1)) | position;
				unsigned int lower = upper + span;
				unsigned int twiddle = position << (LOG2_NUM_OF_SAMPLES - 1 - stage);

				// upper and lower are in different banks; find which is which
				bool upper_in_bank1 = bank_of(upper);
				unsigned int address0 = (upper_in_bank1 ? lower : upper) / 2;
				unsigned int address1 = (upper_in_bank1 ? upper : lower) / 2;

				DTYPE read0_real = bank0_real[address0], read0_imag = bank0_imag[address0];
				DTYPE read1_real = bank1_real[address1], read1_imag = bank1_imag[address1];
				DTYPE upper_real = upper_in_bank1 ? read1_real : read0_real;
				DTYPE upper_imag = upper_in_bank1 ? read1_imag : read0_imag;
				DTYPE lower_real = upper_in_bank1 ? read0_real : read1_real;
				DTYPE lower_imag = upper_in_bank1 ? read0_imag : read1_imag;

				// lower * W
				DTYPE w_real = twiddle_real[twiddle];
				DTYPE w_imag = twiddle_imag[twiddle];
				DTYPE product_real = lower_real * w_real - lower_imag * w_imag;
				DTYPE product_imag = lower_real * w_imag + lower_imag * w_real;

				DTYPE new_upper_real = upper_real + product_real;
				DTYPE new_upper_imag = upper_imag + product_imag;
				DTYPE new_lower_real = upper_real - product_real;
				DTYPE new_lower_imag = upper_imag - product_imag;

				bank0_real[address0] = upper_in_bank1 ? new_lower_real : new_upper_real;
				bank0_imag[address0] = upper_in_bank1 ? new_lower_imag : new_upper_imag;
				bank1_real[address1] = upper_in_bank1 ? new_upper_real : new_lower_real;
				bank1_imag[address1] = upper_in_bank1 ? new_upper_imag : new_lower_imag;
			}
		}

	store_loop:
		for (unsigned int sample = 0; sample < NUM_OF_SAMPLES; sample++) {
#pragma HLS pipeline II=1
			bool in_bank1 = bank_of(sample);
			complex_sample value;
			value.real = in_bank1 ? bank1_real[sample / 2] : bank0_real[sample / 2];
			DTYPE imag = in_bank1 ? bank1_imag[sample / 2] : bank0_imag[sample / 2];
			value.imag = inverse ? -imag : imag;
			output.write(value);
		}
	}
}
