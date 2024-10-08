// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#include <gpc/gpc.h>
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

// Calculates inverse DFT using synthesis
// static void dft_inverse(double xs[], const double dft_re[], const double dft_im[], size_t n)
// {
//     memset(xs, 0, n * sizeof xs[0]);
//
//     for (size_t i = 0; i < n; ++i) {
//         for (size_t j = 1; j <= n/2; ++j) {
//             xs[i] += dft_re[i]*cos(2*FFT_PI*i*j/n) - dft_im[i]*sin(2*FFT_PI*i*j/n);
//         }
//     }
// }

// Simple testing filter, window_size is expected to be small
void moving_average(double ys[], const double xs[], size_t length, size_t window_size)
{
    gp_assert(window_size & 1);
    gp_assert(window_size >= 3, "%zu", window_size);

    // Extrapolate beginning
    for (int i = 0; i < (int)window_size/2; ++i) {
        double sum = 0.;
        for (int j = i - (int)window_size/2; j <= i + (int)window_size/2; ++j) {
            if (j < 0)
                sum += xs[0] + j*(xs[1] - xs[0]);
            else
                sum += xs[j];
        }
        ys[i] = sum/window_size;
    }
    // Calculate middle part
    for (size_t i = window_size/2; i < length - window_size/2; ++i) {
        double sum = 0.;
        for (size_t j = i - window_size/2; j <= i + window_size/2; ++j)
            sum += xs[j];
        ys[i] = sum/window_size;
    }
    // Extrapolate end
    for (int i = length - (int)window_size/2; i <= (int)length; ++i) {
        double sum = 0.;
        for (int j = i - (int)window_size/2; j <= i + (int)window_size/2; ++j) {
            if (j < (int)length)
                sum += xs[j];
            else
                sum += xs[length - 1] + (j - length)*(xs[length - 1] - xs[length - 2]);
        }
        ys[i] = sum/window_size;
    }
}

void print_signal(const char* name, double signal[], size_t length)
{
    #define HEIGHT 15 // should be odd
    #define WIDTH 100

    static double* interpolated = NULL;
    static double* filtered     = NULL;

    char graphic[HEIGHT][WIDTH];
    memset(graphic, ' ', sizeof graphic);

    // Draw edges
    memset(graphic[0], '-', sizeof graphic[0]);
    memset(graphic[HEIGHT / 2], '-', sizeof graphic[HEIGHT / 2]); // zero line
    memset(graphic[HEIGHT - 1], '-', sizeof graphic[HEIGHT - 1]);
    for (size_t y = 0; y < HEIGHT; ++y)
        graphic[y][0] = graphic[y][WIDTH - 1] = '|';
    graphic[0][0] = graphic[0][WIDTH - 1] = graphic[HEIGHT - 1][0] = graphic[HEIGHT - 1][WIDTH - 1] = '+';

    if (length < WIDTH) { // interpolate
        size_t interpolated_length = length;
        while (interpolated_length < WIDTH)
            interpolated_length *= 2;
        size_t k = interpolated_length/length;

        interpolated = realloc(interpolated, 2 * interpolated_length * sizeof interpolated[0]);
        for (size_t i = 0; i < length - 1; ++i) {
            for (size_t j = 0; j < k; ++j)
                interpolated[k*i + j] = (1. - (double)j/k)*signal[i] + ((double)j/k)*signal[i + 1];
        }
        // Extrapolate last samples
        double extrapolated = signal[length - 1] + (signal[length - 1] - signal[length - 2]);
        for (size_t j = 0; j < k; ++j)
            interpolated[k*(length - 1) + j] = (1. - (double)j/k)*signal[length - 1] + ((double)j/k)*extrapolated;
        signal = interpolated;
        length = interpolated_length;
    }

    if (length != WIDTH) { // filter noise for decimation
        size_t window_size = 2*length/WIDTH;
        window_size |= (window_size + 1) & 1; // make odd
        filtered = realloc(filtered, length * sizeof filtered[0]);
        moving_average(filtered, signal, length, window_size);
        signal = filtered;
    }
    // Draw graphic
    for (size_t x = 0; x < WIDTH; ++x) {
        double val = length != WIDTH ?
            signal[x*length/WIDTH] // decimate
          : signal[x];

        if (val == 1.)
            graphic[0][x] = '+';
        else if (val > 1.)
            graphic[0][x] = '~';
        else if (val < -1.)
            graphic[HEIGHT - 1][x] = '~';
        else {
            size_t y = HEIGHT*(.5*val + .5);
            graphic[HEIGHT - 1 - y][x] = '+';
        }
    }
    // Render
    puts(name);
    for (size_t y = 0; y < HEIGHT; ++y) {
        for (size_t x = 0; x < WIDTH; ++x) {
            putchar(graphic[y][x]);
        }
        puts("");
    }
    puts("");
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
                gp_assert(res[i] == sorted[i] && ims[i] == sorted[i],
                    "%zu", i, "%g", res[i], "%g", ims[i], "%g", sorted[i]);
        }
        gp_test("FFT"); {
            #define SIZE 128
            double dft_re[SIZE];
            double dft_im[SIZE];
            double fft_re[SIZE];
            double fft_im[SIZE] = {0};
            GPRandomState rs = gp_new_random_state(time(NULL));
            for (size_t i = 0; i < SIZE; ++i)
                fft_re[i] = 2*gp_frandom(&rs) - 1.;

            dft(dft_re, dft_im, fft_re, SIZE);
            fft_complex(fft_re, fft_im, SIZE);

            print_signal("dft_re", dft_re, SIZE);
            print_signal("dft_im", dft_im, SIZE);
            print_signal("fft_re", fft_re, SIZE);
            print_signal("fft_im", fft_im, SIZE);

            for (size_t i = 0; i < SIZE; ++i) {
                gp_assert(gp_fapprox(fft_re[i], dft_re[i], 0.00001), "%zu", i, "%g", fft_re[i], "%g", dft_re[i]);
                gp_assert(gp_fapprox(fft_im[i], dft_im[i], 0.00001), "%zu", i, "%g", fft_im[i], "%g", dft_im[i]);
            }
        }
    } // gp_suite("FFT")
}
