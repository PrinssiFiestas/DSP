// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#ifndef IIR_H_INCLUDED
#define IIR_H_INCLUDED 1

#include <stddef.h>

#if __cplusplus
extern "C" {
#endif

#if __GNUC__
#define IIR_NONNULL_ARGS(...) __attribute__((nonnull (__VA_ARGS__)))
#else
#define IIR_NONNULL_ARGS(...)
#endif


// ----------------------------------------------------------------------------
//
//           API REFERENCE
//
// ----------------------------------------------------------------------------


/** Filter state array sizes. */
#define IIR_POLES(N) ((N/*previous samples*/) + 1/*current sample*/)

/** Filter state and coefficients.
 * This should always be used as an array large enough to hold the state and
 * coefficients for any given type of filter. For every 6 dB/octave, 1 pole is
 * needed, so a 12 dB/octave filter has to be of type `IIRFilter[IIR_POLES(2)]`.
 * Additionally, in case of stereo input, a 2-dimensional array is needed, so
 * 24 dB/octave stereo filter would be of type `IIRFilter[2][IIR_POLES(4)]`.
 *     The 0th element corresponds to the current sample, 1th corresponds to the
 * previous sample, 2th to the one before, etc.
 */
typedef struct iir_filter
{
    double x; /**< input state  */
    double y; /**< output state */
    double a; /**< input coefficient  */
    double b; /**< output coefficient */
} IIRFilter;

/** Non-linearities
 * Like @ref IIRFilter, this should be used as an array with the same size as
 * the filter, so for `IIRFilter[IIR_POLES(N)]` use
 * `IIRNonlinearities[IIR_POLES(N)]`. Similarly, the 0th element corresponds to
 * the current sample, 1th to one before, etc. Use with @ref iir_apply_filter().
 */
typedef struct iir_nonlinearities
{
    double(*fx)(double); /**< applied to inputs  */
    double(*fy)(double); /**< applied to outputs */
} IIRNonlinearities;

/** @defgroup SampleRate Global Sample Rate
 * Most filters use global sample rate with the assumption of them only being
 * called from the real-time audio thread. If filters are used in multiple
 * threads, or smaple rate changes due to oversampling, or multiple sample rates
 * are used for any other reason, use @ref iir_apply_filter() instead. @{
 */
void   iir_set_sample_rate(double); /**< also sets iir_sample_time().*/
double iir_sample_rate(void);
double iir_sample_time(void);       /**< 1./iir_sample_rate() */
/** @} */

/** @defgroup BasicFilters Basic Precise Linear Filters
 * @p freq can be in Hz or normalized range from 0.0 to 1.0. Uses global sample
 * rate when using Hz. If this is not desired, use @ref iir_apply_filter(). @{
 */
double iir_low_pass6(IIRFilter[IIR_POLES(1)], double input, double freq) IIR_NONNULL_ARGS();
double iir_high_pass6(IIRFilter[IIR_POLES(1)], double input, double freq) IIR_NONNULL_ARGS();
double iir_low_pass12(IIRFilter[IIR_POLES(2)], double input, double freq, double q) IIR_NONNULL_ARGS();
double iir_high_pass12(IIRFilter[IIR_POLES(2)], double input, double freq, double q) IIR_NONNULL_ARGS();
double iir_band_pass12(IIRFilter[IIR_POLES(2)], double input, double freq, double q) IIR_NONNULL_ARGS();
double iir_band_stop12(IIRFilter[IIR_POLES(2)], double input, double freq, double q) IIR_NONNULL_ARGS();
double iir_butterworth_low_pass12(IIRFilter[IIR_POLES(2)], double input, double freq) IIR_NONNULL_ARGS();  /**< Maximally flat pass-band. */
double iir_butterworth_high_pass12(IIRFilter[IIR_POLES(2)], double input, double freq) IIR_NONNULL_ARGS(); /**< Maximally flat pass-band. */
double iir_all_pass6(IIRFilter[IIR_POLES(1)], double input, double freq) IIR_NONNULL_ARGS();
double iir_all_pass12(IIRFilter[IIR_POLES(2)], double input, double freq, double q) IIR_NONNULL_ARGS();
double iir_low_shelving(IIRFilter[IIR_POLES(1)], double input, double freq, double gain) IIR_NONNULL_ARGS();
double iir_high_shelving(IIRFilter[IIR_POLES(1)], double input, double freq, double gain) IIR_NONNULL_ARGS();
double iir_peak(IIRFilter[IIR_POLES(2)], double input, double freq, double q, double gain) IIR_NONNULL_ARGS();
double iir_peak_const_q(IIRFilter[IIR_POLES(2)], double input, double freq, double q, double gain) IIR_NONNULL_ARGS();
double iir_linkwitz_riley_low_pass12(IIRFilter[IIR_POLES(2)], double input, double freq) IIR_NONNULL_ARGS();  /**< -6 dB cutoff point. Useful for cross-overs. */
double iir_linkwitz_riley_high_pass12(IIRFilter[IIR_POLES(2)], double input, double freq) IIR_NONNULL_ARGS(); /**< -6 dB cutoff point. Useful for cross-overs. */
/** @} */

/** @defgroup FastFilters Speed Optimized Linear Filters
 * @p freq will not be accurate, and it may have a small effect on @p damping
 * and vice verca. @p freq shuold be in range 0.0 to 1.0. Use these when
 * precision is not important e.g. musical filtering. @{
 */
double iir_fast_low_pass6(IIRFilter[IIR_POLES(1)], double input, double normalized_freq) IIR_NONNULL_ARGS();
double iir_fast_high_pass6(IIRFilter[IIR_POLES(1)], double input, double normalized_freq) IIR_NONNULL_ARGS();
double iir_fast_low_pass12(IIRFilter[IIR_POLES(2)], double input, double normalized_freq, double damping) IIR_NONNULL_ARGS();
double iir_fast_high_pass12(IIRFilter[IIR_POLES(2)], double input, double normalized_freq, double damping) IIR_NONNULL_ARGS();
/** @} */

/** @defgroup Chebyshevs Chebyshev Linear Filters
 * Chebyshev filtering with even number of @p poles from 2 to 28. @p filter
 * obviously has to be of size `IIR_POLES(@p poles)`.
 *     The filters get more unstable close to 0 Hz and Nyquist with increased
 * poles. It is the responsability of the user to limit @p freq to a stable
 * range. @p freq can be in Hz or from 0.0 to 1.0.
 *     @p ripple should be in range 0.0 to 1.0, where 0.0 corresponds to 0%
 * ripple (Butterworth) and 1.0 corresponds to about 30% ripple. Higher
 * @p ripple gives faster roll-off, but more ripple (duh). @{
 */
double iir_chebyshev_low_pass(IIRFilter[], size_t poles, double input, double freq, double ripple) IIR_NONNULL_ARGS();
double iir_chebyshev_high_pass(IIRFilter[], size_t poles, double input, double freq, double ripple) IIR_NONNULL_ARGS();
/** @} */

// ----------------------------------------------------------------------------
// Generic Filtering

/** Filter with precalculated coefficients and optional nonlinearities.
 * @p filter obviously has to be of size `IIR_POLES(@p poles)`. Same with
 * non-linearities if passed, `NULL` for linear filtering. @p filter must have
 * it's coefficients precalculated. Use a function from @ref Coefficients to
 * calculate coefficients before calling this function.
 *     This library calculates filter output with the recursion equation in form
 * `y0 = a0*x0 + a1*x1 + ... + an*xn + b1*y1 + b2*y2 + ... + bn*yn` as opposed
 * to the difference equation in form of
 * `y0 = a0*x0 + a1*x1 + ... + an*xn - b1*y1 - b2*y2 - ... - bn*yn`. Be mindful
 * of the signs of the b coefficients when writing custom IIR filters.
 */
IIR_NONNULL_ARGS(1)
double iir_apply_filter(
    IIRFilter[],
    size_t poles, // 1 for 6dB/oct, 2 for 12dB/oct, etc.
    double input,
    const IIRNonlinearities optional[]);

/** @defgroup Coefficients Filter Coefficient Precalculation
 * Use these to define filter behaviour when passing filters to
 * @ref iir_apply_filter(). @p sample_time should be the reciprocal of the
 * desired sample rate. @{
 */
void iir_coeffs_low_pass6(IIRFilter[IIR_POLES(1)], double sample_time, double freq) IIR_NONNULL_ARGS();
void iir_coeffs_high_pass6(IIRFilter[IIR_POLES(1)], double sample_time, double freq) IIR_NONNULL_ARGS();
void iir_coeffs_low_pass12(IIRFilter[IIR_POLES(2)], double sample_time, double freq, double q) IIR_NONNULL_ARGS();
void iir_coeffs_high_pass12(IIRFilter[IIR_POLES(2)], double sample_time, double freq, double q) IIR_NONNULL_ARGS();
void iir_coeffs_band_pass12(IIRFilter[IIR_POLES(2)], double sample_time, double freq, double q) IIR_NONNULL_ARGS();
void iir_coeffs_band_stop12(IIRFilter[IIR_POLES(2)], double sample_time, double freq, double q) IIR_NONNULL_ARGS();
void iir_coeffs_butterworth_low_pass12(IIRFilter[IIR_POLES(2)], double sample_time, double freq) IIR_NONNULL_ARGS();  /**< Maximally flat pass-band. */
void iir_coeffs_butterworth_high_pass12(IIRFilter[IIR_POLES(2)], double sample_time, double freq) IIR_NONNULL_ARGS(); /**< Maximally flat pass-band. */
void iir_coeffs_all_pass6(IIRFilter[IIR_POLES(1)], double sample_time, double freq) IIR_NONNULL_ARGS();
void iir_coeffs_all_pass12(IIRFilter[IIR_POLES(2)], double sample_time, double freq, double q) IIR_NONNULL_ARGS();
void iir_coeffs_low_shelving(IIRFilter[IIR_POLES(1)], double sample_time, double freq, double gain) IIR_NONNULL_ARGS();
void iir_coeffs_high_shelving(IIRFilter[IIR_POLES(1)], double sample_time, double freq, double gain) IIR_NONNULL_ARGS();
void iir_coeffs_peak(IIRFilter[IIR_POLES(2)], double sample_time, double freq, double q, double gain) IIR_NONNULL_ARGS();
void iir_coeffs_peak_const_q(IIRFilter[IIR_POLES(2)], double sample_time, double freq, double q, double gain) IIR_NONNULL_ARGS();
void iir_coeffs_linkwitz_riley_low_pass12(IIRFilter[IIR_POLES(2)], double sample_time, double freq) IIR_NONNULL_ARGS();  /**< -6 dB cutoff point. Useful for cross-overs. */
void iir_coeffs_linkwitz_riley_high_pass12(IIRFilter[IIR_POLES(2)], double sample_time, double freq) IIR_NONNULL_ARGS(); /**< -6 dB cutoff point. Useful for cross-overs. */

void iir_coeffs_fast_low_pass6(IIRFilter[IIR_POLES(1)], double normalized_freq) IIR_NONNULL_ARGS();
void iir_coeffs_fast_high_pass6(IIRFilter[IIR_POLES(1)], double normalized_freq) IIR_NONNULL_ARGS();
void iir_coeffs_fast_low_pass12(IIRFilter[IIR_POLES(2)], double normalized_freq, double damping) IIR_NONNULL_ARGS();
void iir_coeffs_fast_high_pass12(IIRFilter[IIR_POLES(2)], double normalized_freq, double damping) IIR_NONNULL_ARGS();

void iir_coeffs_chebyshev_low_pass(IIRFilter[], size_t poles, double sample_time, double freq, double ripple) IIR_NONNULL_ARGS();
void iir_coeffs_chebyshev_high_pass(IIRFilter[], size_t poles, double sample_time, double freq, double ripple) IIR_NONNULL_ARGS();
/** @} */


// ----------------------------------------------------------------------------
//
//          END OF API REFERENCE
//
//          Code below is for internal usage and may change without notice.
//
// ----------------------------------------------------------------------------


#if __cplusplus
} // extern "C"
#endif

#endif // IIR_H_INCLUDED
