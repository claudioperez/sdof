#include <math.h>

#if defined(_WIN32)
#include <Python.h>
PyMODINIT_FUNC PyInit__fsdof(void) {}
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif


struct generalized_alpha {
  double alpha_m,
         alpha_f,
         beta,
         gamma;
} CONF = {1.0, 1.0, 0.25, 0.5};

struct SDOF_Peaks {
    double max_displ,
           max_accel,
           time_max_accel;
};


EXPORT int
fsdof_peaks(struct generalized_alpha* conf,
    double M, double C, double K,
    double scale, int n, double *p, double dt,
    struct SDOF_Peaks *response)
{ 
    conf = &CONF;
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

    double time   = 0.0;
    double   ua,
             va,
             aa,
             u[2],
             v[2],
             a[2];


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

      //  UPDATE(struct *model model, double du)
      u[pres] += du;
      v[pres] += c2*du;
      a[pres] += c3*du;
      
      // COMMIT
      if (fabs(u[pres]) > response->max_displ) {
          response->max_displ = fabs(u[pres]);
      }
      if (fabs(a[pres]) > response->max_accel) {
          response->max_accel = fabs(a[pres]);
          response->time_max_accel = i*dt;
      }
    }
    return 1;
}


EXPORT int
fsdof_integrate(struct generalized_alpha* conf,
    double M, double C, double K,
    double scale, int n, double *p, double dt,
    double *response)
{ 
    conf = &CONF;
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

    double time   = 0.0;
    double  ua,
            va,
            aa,
            *u = &response[0],
            *v = &response[n],
            *a = &response[2*n];


    int i = 0,
        past = -1,
        pres =  0;

//  u[pres] = 0.0;
//  v[pres] = 0.0;
    a[pres] = (p[i] - C*v[pres] - K*u[pres])/M;

    for (i = 1; i < n; i++) {
      ++past;
      ++pres;

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

