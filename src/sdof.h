//
/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2022-2023, Claudio M. Perez
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef SDOF_H
#define SDOF_H

#if defined(_WIN32)
#  define SDOF_EXPORT __declspec(dllexport)

#elif defined(__EMSCRIPTEN__)
#  include <stdlib.h>
#  include <emscripten.h>
#  define SDOF_EXPORT EMSCRIPTEN_KEEPALIVE

#else // *NIXs
#  define SDOF_EXPORT
#endif

#ifdef NEW_API

  /* struct peak_value */
  /* struct  */
  struct sdof_elasticity {
  };

  sdof_peak(int n, double* p, double dt, sdof_elasticity, struct sdof_conf*);
  sdof_hist(int n, double* p, double dt, sdof_elasticity, struct sdof_conf*);
  sdof_accel_hist();
  sdof_accel_peak();
  sdof_veloc_hist();
  sdof_veloc_peak();
  sdof_displ_hist();
  sdof_displ_peak();

#endif

// Parameters for the generalized alpha method.
struct sdof_alpha {
  double alpha_m,
         alpha_f,
         beta,
         gamma;
};

/**
 * Struct to return peak response quantities by value.
 *
 * Members:
 *     max_displ: Peak displacement
 *     max_veloc: Peak velocity
 *     max_accel: Peak acceleration
 *
 */
struct sdof_peaks {
    double max_displ,
           max_veloc,
           max_accel;
};


SDOF_EXPORT struct sdof_peaks
sdof_integrate_peaks_2(struct sdof_alpha* conf,
    double M, double C, double K,
    double scale, int n, const double *p, double dt);

#endif // SDOF_H
