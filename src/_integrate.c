/*
 * Copyright (c) 2022-2023 Claudio Perez
 */

#include "sdof.h"
#include <math.h>

#if defined(_WIN32)
   // Python headers are only required if compiling
   // into a dynamic library that will be loaded by
   // Python on Windows
#  include <Python.h>
   PyMODINIT_FUNC PyInit__integrate(void) {}
#  define EXPORT __declspec(dllexport)

#elif defined(__EMSCRIPTEN__)
#  include <stdlib.h>
#  include <emscripten.h>
#  define EXPORT EMSCRIPTEN_KEEPALIVE

#else // *NIXs
#  define EXPORT
#endif

// Default parameters
struct sdof_alpha CONF = {1.0, 1.0, 0.25, 0.5};


// Main linear integrator. Same as sdof_integrate_0, but operates on transposed
// data. This is better with the cache and faster.
EXPORT int
sdof_integrate(struct sdof_alpha* conf,
    double M, double C, double K,
    double scale, int n, double *p, double dt,
    double *response)
{ 
    const double gamma   = conf->gamma;
    const double beta    = conf->beta;
    const double alpha_m = conf->alpha_m;
    const double alpha_f = conf->alpha_f;

    const double c1 = 1.0;
    const double c2 = gamma/(beta*dt);
    const double c3 = 1.0/(beta*dt*dt);

    const double a1 =     (1.0 -     gamma/beta);
    const double a2 =  dt*(1.0 - 0.5*gamma/beta);
    const double a3 = -1.0/(beta*dt);
    const double a4 =  1.0 - 0.5/beta;

    const double ki = alpha_f*c1*K + alpha_f*c2*C + alpha_m*c3*M;

    double  va,
            aa,
            *u = &response[0],
            *v = &response[1],
            *a = &response[2];

    int i = 0;
    const int past = -3,
              pres =  0;

    // NOTE: The first row of the response array
    // is expected to be initialized!
    //  u[pres] = 0.0;
    //  v[pres] = 0.0;
    a[pres] = (p[i] - C*v[pres] - K*u[pres])/M;

    for (i = 1; i < n; i++) {
      // Move current state pointers forward
      u += 3; v += 3; a += 3;

      // Predictor
      u[pres] = u[past];
      v[pres] = a1*v[past] + a2*a[past];
      a[pres] = a4*a[past] + a3*v[past];

      va = (1.0 - alpha_f)*v[past] + alpha_f*v[pres];
      aa = (1.0 - alpha_m)*a[past] + alpha_m*a[pres];

      //
      // SOLVE
      //
      double pi = (scale*p[i] - C*va - M*aa - K*u[pres]);
      double du = pi / ki;

      //  
      //  UPDATE
      //  
      u[pres] += du;
      v[pres] += c2*du;
      a[pres] += c3*du;
    }
    return 1;
}


#define PAST -3
#define PRES  0
#define UR_STEP(index) do {\
      u += 3; v += 3; a += 3; \
      u[PRES] = u[PAST]; \
      v[PRES] = a1*v[PAST] + a2*a[PAST];      \
      a[PRES] = a4*a[PAST] + a3*v[PAST];      \
      va = (1.0 - alpha_f)*v[PAST] + alpha_f*v[PRES]; \
      aa = (1.0 - alpha_m)*a[PAST] + alpha_m*a[PRES]; \
      /* SOLVE  */ \
      double pi = (scale*p[(index)] - C*va - M*aa - K*u[PRES]);     \
      double du = pi / ki;      \
      /*  UPDATE */             \
      u[PRES] += du;            \
      v[PRES] += c2*du;         \
      a[PRES] += c3*du;         \
    } while (0);

EXPORT int
sdof_integrate_unrolled(struct sdof_alpha* conf,
    const double M, const double C, const double K,
    double scale, const int n, double *p, double dt,
    double *response)
{ 
    const double gamma   = conf->gamma;
    const double beta    = conf->beta;
    const double alpha_m = conf->alpha_m;
    const double alpha_f = conf->alpha_f;

    const double c1 = 1.0;
    const double c2 = gamma/(beta*dt);
    const double c3 = 1.0/(beta*dt*dt);

    const double a1 =     (1.0 -     gamma/beta);
    const double a2 =  dt*(1.0 - 0.5*gamma/beta);
    const double a3 = -1.0/(beta*dt);
    const double a4 =  1.0 - 0.5/beta;

    const double ki = alpha_f*c1*K + alpha_f*c2*C + alpha_m*c3*M;

    double  va,
            aa,
            *u = &response[0],
            *v = &response[1],
            *a = &response[2];

    int i = 0;

    // NOTE: The first row of the response array
    // is expected to be initialized!
    //  u[pres] = 0.0;
    //  v[pres] = 0.0;
    a[PRES] = (p[i] - C*v[PRES] - K*u[PRES])/M;

#define UR_INCR 8
    for (i = 1; i < n-UR_INCR; i+=UR_INCR) {
      UR_STEP(0+i);
      UR_STEP(1+i);
      UR_STEP(2+i);
      UR_STEP(3+i);
      UR_STEP(4+i);
      UR_STEP(5+i);
      UR_STEP(6+i);
      UR_STEP(7+i);
    }
    for (; i < n; i++)
      UR_STEP(0+i);

    return 1;
}


// Integrate a linear system, tracking only peak values.
EXPORT int
sdof_integrate_peaks(struct sdof_alpha* conf,
    double M, double C, double K,
    double scale, int n, double *p, double dt,
    struct sdof_peaks *response)
{ 
    const double gamma   = conf->gamma;
    const double beta    = conf->beta;
    const double alpha_m = conf->alpha_m;
    const double alpha_f = conf->alpha_f;

    const double c1 = 1.0;
    const double c2 = gamma/(beta*dt);
    const double c3 = 1.0/(beta*dt*dt);

    const double a1 =     (1.0 -     gamma/beta);
    const double a2 =  dt*(1.0 - 0.5*gamma/beta);
    const double a3 = -1.0/(beta*dt);
    const double a4 =  1.0 - 0.5/beta;

    const double ki = alpha_f*c1*K + alpha_f*c2*C + alpha_m*c3*M;

    // double time   = 0.0;
    double       va,   aa,
           u[2], v[2], a[2];

    int i = 0, past = 1, pres = 0;

    u[pres] = 0.0;
    v[pres] = 0.0;
    a[pres] = (p[i] - C*v[pres] - K*u[pres])/M;

    for (i = 1; i < n; i++) {
      past = !past;
      pres = !pres;

      u[pres] = u[past];
      v[pres] = a1*v[past] + a2*a[past];
      a[pres] = a4*a[past] + a3*v[past];

      va = (1-alpha_f)*v[past] + alpha_f*v[pres];
      aa = (1-alpha_m)*a[past] + alpha_m*a[pres];      

      // SOLVE
      double pi = (scale*p[i] - C*va - M*aa - K*u[pres]);
      double du = pi / ki;

      u[pres] += du;
      v[pres] += c2*du;
      a[pres] += c3*du;
      
      // COMMIT
      if (fabs(u[pres]) > response->max_displ) {
          response->max_displ = fabs(u[pres]);
      }
      if (fabs(v[pres]) > response->max_veloc) {
          response->max_veloc = fabs(v[pres]);
      }

      double ar = fabs(a[pres] - p[i]/M);
      if (ar > response->max_accel) {
          response->max_accel = ar;
      }
    }
    return 1;
}

// Alternative implementation to return an
// sdof_peaks struct by value. This is
// used in threading.
EXPORT struct sdof_peaks
sdof_integrate_peaks_2(struct sdof_alpha* conf,
    double M, double C, double K,
    double scale, int n, const double *p, double dt)
{ 
    struct sdof_peaks response = {
              .max_displ      = 0.0,
              .max_veloc      = 0.0,
              .max_accel      = 0.0,
    };
    const double gamma   = conf->gamma;
    const double beta    = conf->beta;
    const double alpha_m = conf->alpha_m;
    const double alpha_f = conf->alpha_f;

    const double c1 = 1.0;
    const double c2 = gamma/(beta*dt);
    const double c3 = 1.0/(beta*dt*dt);

    const double a1 =     (1.0 -     gamma/beta);
    const double a2 =  dt*(1.0 - 0.5*gamma/beta);
    const double a3 = -1.0/(beta*dt);
    const double a4 =  1.0 - 0.5/beta;

    const double ki = alpha_f*c1*K + alpha_f*c2*C + alpha_m*c3*M;

    double       va,   aa,
           u[2], v[2], a[2];

    int i = 0, past = 1, pres = 0;

    u[pres] = 0.0;
    v[pres] = 0.0;
    a[pres] = (p[i] - C*v[pres] - K*u[pres])/M;

    for (i = 1; i < n; i++) {
      past = !past;
      pres = !pres;

      u[pres] = u[past];
      v[pres] = a1*v[past] + a2*a[past];
      a[pres] = a4*a[past] + a3*v[past];

      va = (1-alpha_f)*v[past] + alpha_f*v[pres];
      aa = (1-alpha_m)*a[past] + alpha_m*a[pres];      

      // SOLVE
      double pi = (scale*p[i] - C*va - M*aa - K*u[pres]);
      double du = pi / ki;

      u[pres] += du;
      v[pres] += c2*du;
      a[pres] += c3*du;
      
      // COMMIT
      if (fabs(u[pres]) > response.max_displ) {
          response.max_displ = fabs(u[pres]);
      }
      if (fabs(v[pres]) > response.max_veloc) {
          response.max_veloc = fabs(v[pres]);
      }

      double ar = fabs(a[pres] - p[i]/M);
      if (ar > response.max_accel) {
          response.max_accel = ar;
      }
    }
    return response;
}

// RETAINED FOR EDUCATIONAL PURPOSES ONLY
// This implementation uses a slightly different memory layout
// than sdof_integrate. As expected, this implementation is not
// as cache friendly, and is generally inferior
EXPORT int
sdof_integrate_0(struct sdof_alpha* conf,
    double M, double C, double K,
    double scale, int n, double *p, double dt,
    double *response)
{
    const double gamma   = conf->gamma;
    const double beta    = conf->beta;
    const double alpha_m = conf->alpha_m;
    const double alpha_f = conf->alpha_f;

    const double c1 = 1.0;
    const double c2 = gamma/(beta*dt);
    const double c3 = 1.0/(beta*dt*dt);

    const double a1 =     (1.0 -     gamma/beta);
    const double a2 =  dt*(1.0 - 0.5*gamma/beta);
    const double a3 = -1.0/(beta*dt);
    const double a4 =  1.0 - 0.5/beta;

    const double ki = alpha_f*c1*K + alpha_f*c2*C + alpha_m*c3*M;

    double  va,
            aa,
            *u = &response[0],
            *v = &response[n],
            *a = &response[2*n];


    int i = 0;
    const int    past = -1,
                 pres =  0;

//  u[pres] = 0.0;
//  v[pres] = 0.0;
    a[pres] = (p[i] - C*v[pres] - K*u[pres])/M;

    for (i = 1; i < n; i++) {
      // Move current state pointers forward
      ++u; ++v; ++a;

      u[pres] = u[past];
      v[pres] = a1*v[past] + a2*a[past];
      a[pres] = a4*a[past] + a3*v[past];

      va = (1.0 - alpha_f)*v[past] + alpha_f*v[pres];
      aa = (1.0 - alpha_m)*a[past] + alpha_m*a[pres];


      //
      // SOLVE
      //
//    time += alpha_f*dt;
      double pi = (scale*p[i] - C*va - M*aa - K*u[pres]);
      double du = pi / ki;

      //  
      //  UPDATE(struct *model model, double du)
      //  
      u[pres] += du;
      v[pres] += c2*du;
      a[pres] += c3*du;

      // ua = (1-alpha_f) * u[past] + alpha_f * u[pres];
      // va = (1-alpha_f) * v[past] + alpha_f * v[pres];
      // aa = (1-alpha_m) * a[past] + alpha_m * a[pres];
      
      // 
      // COMMIT
      //
//    time += (1.0-alpha_f)*dt;
    }
    return 1;
}


