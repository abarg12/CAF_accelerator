#ifndef CAF_H_
#define CAF_H_

typedef float DTYPE;

#define LOG2_NUM_OF_SAMPLES 8                    /* Number of FFT stages */
#define NUM_OF_SAMPLES (1 << LOG2_NUM_OF_SAMPLES) /* Samples per signal = FFT size */
#define NUM_OF_FREQ_BINS 256                     /* Doppler bins */

/*
 * Top level function: Cross Ambiguity Function.
 *
 * Inputs are the transmit and receive signals split into real and imaginary
 * parts, the Doppler bin spacing in Hz, and the sample rate in Hz.
 *
 * The output is the range-Doppler map, stored row by row with
 * NUM_OF_SAMPLES rows (delay) and NUM_OF_FREQ_BINS columns (Doppler). It
 * matches calculate_cross_ambiguity() in baseline/caf_baseline.py.
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
 * In place radix-2 FFT. Set inverse to true for the (unscaled) inverse FFT.
 */
void fft(DTYPE X_R[NUM_OF_SAMPLES], DTYPE X_I[NUM_OF_SAMPLES], bool inverse);

#endif
