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
#include "sdof.h"
#include <math.h>

#if defined(_WIN32)
#  include <Python.h>
   PyMODINIT_FUNC PyInit__spectrum(void) {}
#  define EXPORT __declspec(dllexport)
#  define C11THREADS

#elif defined(__EMSCRIPTEN__)
#  include <stdlib.h>
#  include <emscripten.h>
#  define EXPORT EMSCRIPTEN_KEEPALIVE

#else /* *NIXs */
#  define EXPORT
#endif

#ifndef M_PI
# define M_PI (3.14159265358979323846)
#endif
/* Pre-define 2*pi */
#define PI_PI (9.869604401089358)

#if !defined(C11THREADS)
# include <pthread.h>
#else
 #ifdef _WIN32
   // MSVC does not implement C11 threads
   #include "tinycthread.h"
 #else
   #include <threads.h>
 #endif
# define pthread_create(a,b,c,d) thrd_create(a, c, d)
# define pthread_t               thrd_t
# define pthread_exit(a)         thrd_exit(a)
# define pthread_join(a,b)       thrd_join(a,b)
#endif

struct thread_data {
  struct sdof_peaks *response;
  int thread_index;
  double t_slp,
         t_min;

  int stride;
  int count;

  struct sdof_alpha* conf;
  double damp;
  int n; 
  const double *p;
  double dt;
};

static void *
run_peaks(void *thread_data) {
  struct thread_data td = *((struct thread_data*) thread_data);
  struct sdof_peaks *const resp = td.response;
  const int i = td.thread_index;
  const int stride   = td.stride;
  const int count    = td.count;

  const double scale = 1.0,
               mass  = 1.0,
               damp  = td.damp;

  for (int offset = 0 ; offset < count; offset++) {
    double period = td.t_min + td.t_slp*((double)(i*stride+offset));
    double K = 4.0*PI_PI*mass/(period*period);
    double C = 4.0*damp*M_PI/period;
    resp[i*stride + offset] = sdof_integrate_peaks_2(td.conf,   mass, C, K,
                                     scale, td.n, td.p, td.dt);
  }

  pthread_exit(EXIT_SUCCESS);
  return NULL;
}

// Threaded response spectrum
//
EXPORT int
sdof_spectrum(struct sdof_alpha* conf,
              const double* load, const int n, const double dt, 
              const double t_min, const double t_max, const int n_periods,
              const double damp,
              int n_threads,
              struct sdof_peaks *response)
{
  pthread_t *threads = malloc(sizeof(pthread_t)*n_threads);
  struct thread_data *wkspace = malloc(sizeof(struct thread_data)*n_threads);

  double slope = (t_max - t_min)/((double)n_periods);

  int i;
  for (i = 0; i < n_threads-1; i++) {
    wkspace[i].damp   =  damp;
    wkspace[i].n      =     n;
    wkspace[i].dt     =    dt;
    wkspace[i].p      =  load;
    wkspace[i].conf   =  conf;
    wkspace[i].stride = n_periods/n_threads;
    wkspace[i].count  = n_periods/n_threads;

    wkspace[i].t_min  = t_min;
    wkspace[i].t_slp  = slope;
    wkspace[i].response = response;
    wkspace[i].thread_index = i;

    pthread_create(&threads[i], NULL, &run_peaks, (void *)&wkspace[i]);
  }

  wkspace[i].damp   =  damp;
  wkspace[i].n      =     n;
  wkspace[i].dt     =    dt;
  wkspace[i].p      =  load;
  wkspace[i].conf   =  conf;
  wkspace[i].stride = n_periods/n_threads;
  wkspace[i].count  = n_periods/n_threads + n_periods%n_threads;

  wkspace[i].t_min  = t_min;
  wkspace[i].t_slp  = slope;
  wkspace[i].response = response;
  wkspace[i].thread_index = i;

  pthread_create(&threads[i], NULL, &run_peaks, (void *)&wkspace[i]);

  for(int i = 0; i < n_threads; i++)
    pthread_join(threads[i], NULL);

  free(threads);
  free(wkspace);

  return 0;
}
