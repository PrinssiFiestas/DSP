// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#include <gpc/assert.h>
#include "../src/fft.c"
#include <time.h>

// Calculates DFT using correlation
static void dft(double dft_re[], double dft_im[], const double xs[], size_t n)
{
    memset(dft_re, 0, n * sizeof dft_re[0]);
    memset(dft_im, 0, n * sizeof dft_im[0]);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            dft_re[i] += xs[i]*cos(2*FFT_PI*i*j/n);
            dft_im[i] -= xs[i]*sin(2*FFT_PI*i*j/n);
        }
    }
}

int main(void)
{
    gp_suite("FFT"); {
        gp_test("Bit reversal sort"); {
            double res[8]    = { 0, 1, 2, 3, 4, 5, 6, 7 };
            double ims[8]    = { 0, 1, 2, 3, 4, 5, 6, 7 };
            double sorted[8] = { 0, 4, 2, 6, 1, 5, 3, 7 }; // bit patterns reversed
            fft_bit_reversal_sort(res, ims, sizeof res/sizeof res[0]);
            for (size_t i = 0; i < sizeof res/sizeof res[0]; ++i)
                gp_assert(res[i] == sorted[i] && ims[i] == sorted[i], i, res[i], ims[i], sorted[i]);
        }
        gp_test("FFT"); {
            #define SIZE 1024
            double dft_re[SIZE];
            double dft_im[SIZE];
            double fft_re[SIZE];
            double fft_im[SIZE] = {0};
            GPRandomState rs = gp_new_random_state(time(NULL));
            for (size_t i = 0; i < SIZE; ++i)
                fft_re[i] = 2*gp_frandom(&rs) - 1.;

            dft(dft_re, dft_im, fft_re, SIZE);
            fft_complex(fft_re, fft_im, SIZE);

            for (size_t i = 0; i < SIZE; ++i) {
                gp_assert(gp_fapprox(fft_re[i], dft_re[i], 0.00001), i, fft_re[i], dft_re[i]);
                gp_assert(gp_fapprox(fft_im[i], dft_im[i], 0.00001), i, fft_im[i], dft_im[i]);
            }
        }
    } // gp_suite("FFT")
}
