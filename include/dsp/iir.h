// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#ifndef DSP_IIR_INCLUDED
#define DSP_IIR_INCLUDED 1

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

typedef struct iir_coeffs
{
    double a[IIR_MAX_POLES + 1];
    double b[IIR_MAX_POLES + 1];
} IIRCoeffs;

typedef struct iir_filter
{
    IIRCoeffs coeff;
    double x[IIR_MAX_POLES + 1];
    double y[IIR_MAX_POLES + 1];
} IIRFilter;

typedef struct iir_nonlinearities
{
    double(*fx[IIR_MAX_POLES + 1])(double);
    double(*fy[IIR_MAX_POLES + 1])(double);
} IIRNonlinearities;

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

double iir_fast_lp6(IIRFilter*, double input, double normalized_freq) IIR_NONNULL_ARGS();
double iir_fast_hp6(IIRFilter*, double input, double normalized_freq) IIR_NONNULL_ARGS();
double iir_fast_lp12(IIRFilter*, double input, double normalized_freq, double damping) IIR_NONNULL_ARGS();
double iir_fast_hp12(IIRFilter*, double input, double normalized_freq, double damping) IIR_NONNULL_ARGS();

/** Filter with precalculated coefficients and optional nonlinearities.*/
double iir_apply_filter(IIRFilter*, double input, IIRNonlinearities* optional) IIR_NONNULL_ARGS(1);

/** Create table of nonlinearities with all elements being @p f.*/
IIRNonlinearities iir_nonlinearities(double(*f)(double x)) IIR_NONNULL_ARGS();

IIRCoeffs iir_coeffs_low_pass6(double sample_rate, double freq);
IIRCoeffs iir_coeffs_high_pass6(double sample_rate, double freq);
IIRCoeffs iir_coeffs_low_pass12(double sample_rate, double freq, double q);
IIRCoeffs iir_coeffs_high_pass12(double sample_rate, double freq, double q);
IIRCoeffs iir_coeffs_band_pass12(double sample_rate, double freq, double q);
IIRCoeffs iir_coeffs_band_stop12(double sample_rate, double freq, double q);
IIRCoeffs iir_coeffs_butterworth_low_pass12(double sample_rate, double freq);
IIRCoeffs iir_coeffs_butterworth_high_pass12(double sample_rate, double freq);
IIRCoeffs iir_coeffs_all_pass6(double sample_rate, double freq);
IIRCoeffs iir_coeffs_all_pass12(double sample_rate, double freq, double q);
IIRCoeffs iir_coeffs_low_shelving(double sample_rate, double freq, double gain_db);
IIRCoeffs iir_coeffs_high_shelving(double sample_rate, double freq, double gain_db);
IIRCoeffs iir_coeffs_peak(double sample_rate, double freq, double q, double gain_db);
IIRCoeffs iir_coeffs_peak_const_q(double sample_rate, double freq, double q, double gain_db);
IIRCoeffs iir_coeffs_linkwitz_riley_lp12(double sample_rate, double freq);
IIRCoeffs iir_coeffs_linkwitz_riley_hp12(double sample_rate, double freq);

IIRCoeffs iir_coeffs_fast_lp6(double sample_rate, double normalized_freq);
IIRCoeffs iir_coeffs_fast_hp6(double sample_rate, double normalized_freq);
IIRCoeffs iir_coeffs_fast_lp12(double sample_rate, double normalized_freq, double damping);
IIRCoeffs iir_coeffs_fast_hp12(double sample_rate, double normalized_freq, double damping);


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

#endif // DSP_IIR_INCLUDED
