//

/*
 * BSD 2-Clause License
 *
 * Copyright (c) 2022-2024, Claudio M. Perez
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
#include <stddef.h>

#if defined(_WIN32)
   // Python header is only required if compiling
   // into a dynamic library that will be loaded by
   // Python on Windows. All we really need is the
   // declaration of PyMODINIT_FUNC
#  include <Python.h>
   PyMODINIT_FUNC PyInit__integrate(void) {}
#endif

// Default parameters
struct sdof_alpha CONF = {1.0, 1.0, 0.25, 0.5};

/**
 * Main linear integrator.
 *
 *
 * Parameters:
 *     conf (struct sdof_alpha*): Struct holding integration parameters.
 *
 *     const double[n] load: Pointer to excitation series
 *     const int n: size of excitation series
 *     const double dt: Time step of excitation data.
 *
 *     M (double): Mass.
 *     C (double): Damping.
 *     K (double): Stiffness
 *
 *     response (double[n][3]): Pointer to the beginning of an n x 3 array of doubles.
 *        displacement at ``i``: ``response[i][0]``
 *        velocity at ``i``:     ``response[i][1]``
 *        acceleration at ``i``: ``response[i][2]``
 *
 * .. note:: The first row of the response array
 *     is expected to be initialized, ie
 *     ::
 *         response[0][0] = u0;
 *         response[0][1] = v0;
 */
SDOF_EXPORT int
sdof_integrate(struct sdof_alpha* conf,
    double M, double C, double K,
    double scale, int n, double *p, double dt,
    double *response)
{ 
    if (conf == NULL)
      conf = &CONF;

    const double gamma   = conf->gamma;
    const double beta    = conf->beta;
    const double alpha_m = conf->alpha_m;
    const double alpha_f = conf->alpha_f;

    const double c1 = 1.0;
    const double c2 = gamma/(beta*dt);
    const double c3 = 1.0/(beta*dt*dt);

    const double b1 =     (1.0 -     gamma/beta);
    const double b2 =  dt*(1.0 - 0.5*gamma/beta);
    const double b3 = -1.0/(beta*dt);
    const double b4 =  1.0 - 0.5/beta;

    const double ki = alpha_f*c1*K + alpha_f*c2*C + alpha_m*c3*M;

    double  ua,
            va,
            aa,
            *u = &response[0],
            *v = &response[1],
            *a = &response[2];

    int i = 0;
    const int past = -3,
              pres =  0;

    a[pres] = (p[i] - C*v[pres] - K*u[pres])/M;

    for (i = 1; i < n; i++) {
      // Move current state pointers forward
      u += 3; v += 3; a += 3;

      // SETUP
      u[pres] = u[past];
      v[pres] = b1*v[past] + b2*a[past];
      a[pres] = b4*a[past] + b3*v[past];

      {
        ua = (1.0 - alpha_f)*u[past] + alpha_f*u[pres];
        va = (1.0 - alpha_f)*v[past] + alpha_f*v[pres];
        aa = (1.0 - alpha_m)*a[past] + alpha_m*a[pres];

        //
        // SOLVE
        //
        double ri = scale*p[i] - (M*aa + C*va + K*ua);
        double du = ri / ki;

        //
        // UPDATE
        //
        u[pres] += c1*du;
        v[pres] += c2*du;
        a[pres] += c3*du;
      }
    }
    return 1;
}

/**
 * Elastic-perfectly plastic
 *
 * - Simo, J.C. and Hughes, T.J.R. (2000) Computational inelasticity. Corr. 2.
 *   print. New York Heidelberg Berlin: Springer (Interdisciplinary applied
 *   mathematics Mechanics and materials, 7).
 *
 * Parameters:
 *   Fy (double): Yield force
 *   a (double): Kinematic hardening ratio
 *
 */
SDOF_EXPORT int
sdof_integrate_plastic(struct sdof_alpha* conf,
    double M, double C, double K, double Fy,
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

    const double b1 =     (1.0 -     gamma/beta);
    const double b2 =  dt*(1.0 - 0.5*gamma/beta);
    const double b3 = -1.0/(beta*dt);
    const double b4 =  1.0 - 0.5/beta;

    const double k0 = alpha_f*c2*C + alpha_m*c3*M;

    double  pa = 0.0, // TODO: Find pa for initial u0
            *u = &response[0],
            *v = &response[1],
            *a = &response[2];

    const int past = -3,
              pres =  0;

    // Plasticity
    const int maxIter = 10;
    double b = 0.00,
           tol   = 1e-12;

    struct {double up, Fy, Hkin;} model = {
          .up    = 0.0,
          .Fy    = Fy,
          .Hkin  = K*b/(1.0 - b)
    };


    int i     = 0;
    double Kt = K;
    a[pres]   = (p[i] - C*v[pres] - pa)/M;

    for (i = 1; i < n; i++) {
      // Move current state pointers forward
      u += 3; v += 3; a += 3;

      // Predictor
      u[pres]   = u[past];
      v[pres]   = b1*v[past] + b2*a[past];
      a[pres]   = b4*a[past] + b3*v[past];
      double ua = (1.0 - alpha_f)*u[past] + alpha_f*u[pres];
      double va = (1.0 - alpha_f)*v[past] + alpha_f*v[pres];
      double aa = (1.0 - alpha_m)*a[past] + alpha_m*a[pres];

      double R = scale*p[i] - (M*aa + C*va + pa);
      double R0 = 1.0; //R;
      if (R0 == 0.0)
          R0 = 1.0;

      double du = 0.0;
      for (int iter = 0; iter <= maxIter; iter++) {

        ua += alpha_f*c1*du;
        va += alpha_f*c2*du;
        aa += alpha_m*c3*du;

        // State determination for pa, Kt
        pa = K*(ua - model.up);
        double ftrial = fabs(pa - model.Hkin*model.up) - model.Fy;

        if (ftrial <= 0) {
          // Elastic step
          Kt = K;

        } else {
          // plastic strain increment
          double dup = ftrial/(K + model.Hkin);

          if (pa < 0) {
            pa += dup*K;
            model.up -= dup;
          } else {
            pa -= dup*K;
            model.up += dup;
          }
          Kt = K*model.Hkin/(K + model.Hkin);
        }

        // SOLVE
        R  = scale*p[i] - (M*aa + C*va + pa);
        du = R/(alpha_f*c1*Kt + k0);

        // UPDATE
        u[pres] += c1*du;
        v[pres] += c2*du;
        a[pres] += c3*du;

        if (fabs(R/R0) < tol)
          break;

      }

    }
    return 1;
}



#define PAST -3
#define PRES  0
#define UR_STEP(index) do {\
      u += 3; v += 3; a += 3; \
      u[PRES] = u[PAST]; \
      v[PRES] = b1*v[PAST] + b2*a[PAST];      \
      a[PRES] = b4*a[PAST] + b3*v[PAST];      \
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

SDOF_EXPORT int
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

    const double b1 =     (1.0 -     gamma/beta);
    const double b2 =  dt*(1.0 - 0.5*gamma/beta);
    const double b3 = -1.0/(beta*dt);
    const double b4 =  1.0 - 0.5/beta;

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


/**
 * Integrate a linear system, tracking only peak values.
 */
SDOF_EXPORT int
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

    const double b1 =     (1.0 -     gamma/beta);
    const double b2 =  dt*(1.0 - 0.5*gamma/beta);
    const double b3 = -1.0/(beta*dt);
    const double b4 =  1.0 - 0.5/beta;

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
      v[pres] = b1*v[past] + b2*a[past];
      a[pres] = b4*a[past] + b3*v[past];

      va = (1-alpha_f)*v[past] + alpha_f*v[pres];
      aa = (1-alpha_m)*a[past] + alpha_m*a[pres];      

      // SOLVE
      double pi = (scale*p[i] - C*va - M*aa - K*u[pres]);
      double du = pi / ki;

      u[pres] +=    du;
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

/**
 * Alternative implementation to return an sdof_peaks struct by value. This is
 * used in _spectrum.c .
 */
SDOF_EXPORT struct sdof_peaks
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

    const double b1 =     (1.0 -     gamma/beta);
    const double b2 =  dt*(1.0 - 0.5*gamma/beta);
    const double b3 = -1.0/(beta*dt);
    const double b4 =  1.0 - 0.5/beta;

    const double c1 = 1.0;
    const double c2 = gamma/(beta*dt);
    const double c3 = 1.0/(beta*dt*dt);

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
      v[pres] = b1*v[past] + b2*a[past];
      a[pres] = b4*a[past] + b3*v[past];

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

