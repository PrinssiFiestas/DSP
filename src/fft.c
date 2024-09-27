// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#include <dsp/iir.h>
#include <stdbool.h>
#include <math.h>

#define FFT_E           2.7182818284590452354
#define FFT_LOG2E       1.4426950408889634074
#define FFT_LOG10E      0.43429448190325182765
#define FFT_LN2         0.69314718055994530942
#define FFT_LN10        2.30258509299404568402
#define FFT_PI          3.14159265358979323846
#define FFT_PI_2        1.57079632679489661923
#define FFT_PI_4        0.78539816339744830962
#define FFT_1_PI        0.31830988618379067154
#define FFT_2_PI        0.63661977236758134308
#define FFT_2_SQRTPI    1.12837916709551257390
#define FFT_SQRT2       1.41421356237309504880
#define FFT_SQRT1_2     0.70710678118654752440

// res = reals, ims = imaginarys

static void fft_bit_reversal_sort(double res[], double ims[], size_t length)
{
    size_t j = length/2;
    for (size_t i = 1; i < length - 2; ++i) {
        if (i < j) {
            double tr = res[j];
            double ti = res[j];
            res[j] = res[i];
            ims[j] = ims[i];
            res[i] = tr;
            ims[j] = ti;
        }
        size_t k = length/2;
        while (k <= j) {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}

void fft_complex(double res[], double ims[], size_t length)
{
    size_t stages = FFT_LOG2E*log(length);
    fft_bit_reversal_sort(res, ims, length);

    for (size_t l = 1; l <= stages; ++l) {
        size_t le  = 1 << l;
        size_t le2 = le/2;
        double ur  = 1.;
        double ui  = 0.;
        double sr  =  cos(FFT_PI/le2);
        double si  = -sin(FFT_PI/le2);
        for (size_t j = 0; j < le2; ++j) { // Calculate sub DFT's
            for (size_t i = j; i < length; i += le) { // Calculate butterflys
                size_t ip = i + le2;
                double tr = res[ip]*ur - ims[ip]*ui;
                double ti = res[ip]*ui + ims[ip]*ur;
                res[ip] = res[i] - tr;
                ims[ip] = ims[i] - ti;
                res[i] += tr;
                ims[i] += ti;
            }
            double tr = ur;
            ur = tr*sr - ui*si;
            ui = tr*si + ui*sr;
        }
    }
}
