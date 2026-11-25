#ifndef CAF_H_
#define CAF_H_

#include "hls_stream.h"

typedef float DTYPE;

#define LOG2_NUM_OF_SAMPLES 8                    /* Number of FFT stages */
#define NUM_OF_SAMPLES (1 << LOG2_NUM_OF_SAMPLES) /* Samples per signal = FFT size */
/* Doppler bins. run_hls.tcl overrides this with -DNUM_OF_FREQ_BINS=<n>
 * (even) so that C/RTL co-simulation finishes quickly. */
#ifndef NUM_OF_FREQ_BINS
#define NUM_OF_FREQ_BINS 256
#endif

/* One complex sample, passed between the dataflow processes */
struct complex_sample {
	DTYPE real;
	DTYPE imag;
};

/*
 * Top level function: Cross Ambiguity Function.
 *
 * Same interface and output as the baseline in hls/baseline/caf.h.
 *
 * Inputs are the transmit and receive signals split into real and imaginary
 * parts, the Doppler bin spacing in Hz, and the sample rate in Hz.
 *
 * The output is the range-Doppler map, stored row by row with
 * NUM_OF_SAMPLES rows (delay) and NUM_OF_FREQ_BINS columns (Doppler). It
 * matches the real part of cross_ambiguity() in baseline/caf_baseline.py.
 */
void cross_ambiguity(
    DTYPE transmit_real[NUM_OF_SAMPLES],
    DTYPE transmit_imag[NUM_OF_SAMPLES],
    DTYPE receive_real[NUM_OF_SAMPLES],
    DTYPE receive_imag[NUM_OF_SAMPLES],
    DTYPE freq_step,
    DTYPE sample_rate,
    DTYPE range_doppler_map[NUM_OF_SAMPLES * NUM_OF_FREQ_BINS]);

/*
 * Streaming FFT. Reads num_of_frames frames of NUM_OF_SAMPLES samples from
 * input and writes the FFT of each frame, in natural order, to output. Set
 * inverse to true for the (unscaled) inverse FFT.
 */
void fft_frames(
    hls::stream<complex_sample> &input,
    hls::stream<complex_sample> &output,
    int num_of_frames,
    bool inverse);

#endif
