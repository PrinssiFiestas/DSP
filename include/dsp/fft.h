// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#ifndef FFT_H_INCLUDED
#define FFT_H_INCLUDED 1

#include <stddef.h>

#if __cplusplus
extern "C" {
#endif

#if __GNUC__
#define FFT_NONNULL_ARGS(...) __attribute__((nonnull (__VA_ARGS__)))
#else
#define FFT_NONNULL_ARGS(...)
#endif


// ----------------------------------------------------------------------------
//
//           API REFERENCE
//
// ----------------------------------------------------------------------------


FFT_NONNULL_ARGS()
void fft_complex(double reals[], double imaginarys[], size_t length);


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

#endif // FFT_H_INCLUDED
