// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#ifndef DSP_IIR_INCLUDED
#define DSP_IIR_INCLUDED 1

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
//			PUBLIC API
//
// ----------------------------------------------------------------------------


#ifndef IIR_MAX_POLES
#define IIR_MAX_POLES 4
#endif

#define IIR_POLES(N) ((N) + 1)

typedef struct iir_coeffs
{
    double a[IIR_MAX_POLES + 1];
    double b[IIR_MAX_POLES + 1];
} IIRCoeffs;

typedef struct iir_filter
{
    double x[IIR_MAX_POLES + 1];
    double y[IIR_MAX_POLES + 1];
    IIRCoeffs coeff;
} IIRFilter;

typedef struct iir_nonlinearities
{
    double(*fx[IIR_MAX_POLES + 1])(double);
    double(*fy[IIR_MAX_POLES + 1])(double);
} IIRNonlinearities;

void   iir_set_sample_rate(double); /**< also sets iir_sample_time().*/
double iir_sample_rate(void);
double iir_sample_time(void); /**< 1./iir_sample_rate() */

double iir_low_pass6(IIRFilter*, double input, double freq) IIR_NONNULL_ARGS();
double iir_high_pass6(IIRFilter*, double input, double freq) IIR_NONNULL_ARGS();
double iir_low_pass12(IIRFilter*, double input, double freq, double q) IIR_NONNULL_ARGS();
double iir_high_pass12(IIRFilter*, double input, double freq, double q) IIR_NONNULL_ARGS();
double iir_band_pass12(IIRFilter*, double input, double freq, double q) IIR_NONNULL_ARGS();
double iir_band_stop12(IIRFilter*, double input, double freq, double q) IIR_NONNULL_ARGS();
double iir_butterworth_low_pass12(IIRFilter*, double input, double freq) IIR_NONNULL_ARGS();
double iir_butterworth_high_pass12(IIRFilter*, double input, double freq) IIR_NONNULL_ARGS();
double iir_all_pass6(IIRFilter*, double input, double freq) IIR_NONNULL_ARGS();
double iir_all_pass12(IIRFilter*, double input, double freq, double q) IIR_NONNULL_ARGS();
double iir_low_shelving(IIRFilter*, double input, double freq, double gain_db) IIR_NONNULL_ARGS();
double iir_high_shelving(IIRFilter*, double input, double freq, double gain_db) IIR_NONNULL_ARGS();
double iir_peak(IIRFilter*, double input, double freq, double q, double gain_db) IIR_NONNULL_ARGS();
double iir_peak_const_q(IIRFilter*, double input, double freq, double q, double gain_db) IIR_NONNULL_ARGS();
double iir_linkwitz_riley_lp12(IIRFilter*, double input, double freq) IIR_NONNULL_ARGS();
double iir_linkwitz_riley_hp12(IIRFilter*, double input, double freq) IIR_NONNULL_ARGS();

double iir_fast_low_pass6(IIRFilter*, double input, double normalized_freq) IIR_NONNULL_ARGS();
double iir_fast_high_pass6(IIRFilter*, double input, double normalized_freq) IIR_NONNULL_ARGS();
double iir_fast_low_pass12(IIRFilter*, double input, double normalized_freq, double damping) IIR_NONNULL_ARGS();
double iir_fast_high_pass12(IIRFilter*, double input, double normalized_freq, double damping) IIR_NONNULL_ARGS();

/** Filter with precalculated coefficients and optional nonlinearities.*/
IIR_NONNULL_ARGS(1)
double iir_apply_filter(
    IIRFilter*,
    size_t poles, // 1 for 6dB/oct, 2 for 12dB/oct, etc.
    double input,
    const IIRNonlinearities* optional);

/** Create table of nonlinearities with all elements being @p f.*/
IIRNonlinearities iir_nonlinearities(double(*f)(double x)) IIR_NONNULL_ARGS();

void iir_coeffs_low_pass6(IIRCoeffs*, double sample_time, double freq);
void iir_coeffs_high_pass6(IIRCoeffs*, double sample_time, double freq);
void iir_coeffs_low_pass12(IIRCoeffs*, double sample_time, double freq, double q);
void iir_coeffs_high_pass12(IIRCoeffs*, double sample_time, double freq, double q);
void iir_coeffs_band_pass12(IIRCoeffs*, double sample_time, double freq, double q);
void iir_coeffs_band_stop12(IIRCoeffs*, double sample_time, double freq, double q);
void iir_coeffs_butterworth_low_pass12(IIRCoeffs*, double sample_time, double freq);
void iir_coeffs_butterworth_high_pass12(IIRCoeffs*, double sample_time, double freq);
void iir_coeffs_all_pass6(IIRCoeffs*, double sample_time, double freq);
void iir_coeffs_all_pass12(IIRCoeffs*, double sample_time, double freq, double q);
void iir_coeffs_low_shelving(IIRCoeffs*, double sample_time, double freq, double gain_db);
void iir_coeffs_high_shelving(IIRCoeffs*, double sample_time, double freq, double gain_db);
void iir_coeffs_peak(IIRCoeffs*, double sample_time, double freq, double q, double gain_db);
void iir_coeffs_peak_const_q(IIRCoeffs*, double sample_time, double freq, double q, double gain_db);
void iir_coeffs_linkwitz_riley_low_pass12(IIRCoeffs*, double sample_time, double freq);
void iir_coeffs_linkwitz_riley_high_pass12(IIRCoeffs*, double sample_time, double freq);

void iir_coeffs_fast_low_pass6(IIRCoeffs*, double normalized_freq);
void iir_coeffs_fast_high_pass6(IIRCoeffs*, double normalized_freq);
void iir_coeffs_fast_low_pass12(IIRCoeffs*, double normalized_freq, double damping);
void iir_coeffs_fast_high_pass12(IIRCoeffs*, double normalized_freq, double damping);


// ----------------------------------------------------------------------------
//
//          END OF API REFERENCE
//
//          Code below is for internal usage and may change without notice.
//
// ----------------------------------------------------------------------------

typedef struct filter
{
    double x;
    double y;
    double a;
    double b;
} Filter;

IIR_NONNULL_ARGS()
double low_pass12(Filter filter[3], double input, double freq, double q);

IIR_NONNULL_ARGS()
double fast_low_pass6(Filter filter[2], double input, double freq);

IIR_NONNULL_ARGS()
double fast_low_pass12(Filter filter[3], double input, double freq, double damping);

IIR_NONNULL_ARGS()
double chebyshev_low_pass(Filter filter[], size_t poles, double input, double freq, double ripple);

IIR_NONNULL_ARGS()
double chebyshev_high_pass(Filter filter[], size_t poles, double input, double freq, double ripple);


#if __cplusplus
} // extern "C"
#endif

#endif // DSP_IIR_INCLUDED
