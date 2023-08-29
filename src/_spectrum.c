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
#include "sdof.h"
#include <math.h>
#include <stdlib.h>

#if defined(_WIN32)
#  define C11THREADS
   // Python header is only required if compiling
   // into a dynamic library that will be loaded by
   // Python on Windows. All we really need is the
   // declaration of PyMODINIT_FUNC.
   //
   // Apart from the next two lines, this library has
   // no connection to Python.
#  include <Python.h>
   PyMODINIT_FUNC PyInit__spectrum(void) {}
#endif

#ifndef M_PI
# define M_PI (3.14159265358979323846)
#endif

// Pre-compute 2*pi
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

  // Translate some pthread functions to C11's thrd API
# define pthread_create(a,b,c,d) thrd_create(a, c, d)
# define pthread_t               thrd_t
# define pthread_exit(a)         thrd_exit(a)
# define pthread_join(a,b)       thrd_join(a,b)
#endif

struct sdof_thread {
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
run_peaks(void *sdof_thread) {
  const struct sdof_thread td = *((struct sdof_thread*) sdof_thread);
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

/**
 * Threaded response spectrum over regularly spaced periods.
 *
 * This function spawns `n_threads` threads to perform a total
 * of `n_periods` integrations over regularly spaced periods.
 *
 * Parameters:
 *     struct sdof_alpha* conf: .
 *
 *     const double[n] load: Pointer to excitation series
 *     const int n: size of excitation series
 *     const double dt: Time step of excitation data.
 *
 *     const double t_min: Start period
 *     const double t_max: End period
 *     const int n_periods: Number of periods to integrate.
 *    
 *     const double damp: .
 *     int n_threads: .
 *     struct sdof_peaks *response
 *
 */
SDOF_EXPORT int
sdof_spectrum(struct sdof_alpha* conf,
              const double* load, const int n, const double dt, 
              const double t_min, const double t_max, const int n_periods,
              const double damp,
              int n_threads,
              struct sdof_peaks *response)
{
  pthread_t *threads = malloc(sizeof(pthread_t)*n_threads);
  struct sdof_thread *wkspace = malloc(sizeof(struct sdof_thread)*n_threads);

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
