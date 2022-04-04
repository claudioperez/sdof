#include <stdio.h>
#include <stdlib.h>
#include <math.h>


struct step {double u, v, a;};

struct load {
  int n;
  double dt;
  double p[];
};

struct model {
  double time;
  struct step *past, *pres; 
} model;

struct conf {
  struct model* model;
  double alpha_m, alpha_f, beta, gamma;
} conf = {&model, 1.0, 1.0, 1.0/6.0, 0.5};




int generalized_alpha(struct conf* conf, 
    double M, double C, double K,
    double scale, int n, double p[n], double dt)
{ 
    double gamma   = conf->gamma;
    double beta    = conf->beta;
    double alpha_m = conf->alpha_m;
    double alpha_f = conf->alpha_f;

    double c1 = 1.0;
    double c2 = gamma/(beta*dt);
    double c3 = 1.0/(beta*dt*dt);

    // step coefficients
    double a1 =     (1.0 -     gamma/beta);
    double a2 =  dt*(1.0 - 0.5*gamma/beta);
    double a3 = -1.0/(beta*dt);
    double a4 =  1.0 - 0.5/beta;

    double ki = alpha_f*c1*K + alpha_f*c2*C + alpha_m*c3*M;

    double time = 0.0;
    double   ua   = 0.0,
             va   = 0.0,
             aa   = 0.0,
             u[2] = {0.0, 0.0},
             v[2] = {0.0, 0.0},
             a[2] = {0.0, 0.0};


    int i = 0, past = 1, pres = 0;

    u[pres] = 0.0;
    v[pres] = 0.0;
    a[pres] = (p[i] - C*v[pres] - K*u[pres])/M;
    printf("%lf\t%lf\t%lf\t%lf\n", p[i], u[pres], v[pres], a[pres]);

    for (i = 1; i < n; i++) {
      past = !past;
      pres = !pres;

      u[pres] = u[past];
      v[pres] = a1*v[past] + a2*a[past];
      a[pres] = a4*a[past] + a3*v[past];

      // determine the velocities at t+alpha_f*dt
      // ua = 
      va = (1-alpha_f) * v[past] + alpha_f * v[pres];
      aa = (1-alpha_m) * a[past] + alpha_m * a[pres];
      

      //
      // SOLVE
      //
      time += alpha_f*dt;
      double pi = (scale*p[i] - C*va - M*aa - K*u[pres]);
      double dx = pi / ki;

      /*
       * UPDATE(struct *model model, double dx)
       */
      //  determine the response at t+dt
      u[pres] += dx;
      v[pres] += c2*dx;
      a[pres] += c3*dx;

      // determine displacement and velocity at t+alpha_f*dt
      ua = (1-alpha_f) * u[past] + alpha_f * u[pres];
      va = (1-alpha_f) * v[past] + alpha_f * v[pres];
      aa = (1-alpha_m) * a[past] + alpha_m * a[pres];
      
      // update the response at the DOFs
      // theModel->setResponse(*da[pres],*va,*a[pres]); 

      /* 
       * COMMIT(void)
       */
      printf("%lf\t%lf\t%lf\t%lf\n", pi, u[pres], v[pres], a[pres]);
      time += (1.0-alpha_f)*dt;
    }
    return 1;
}



int read_load(FILE* file, int n, double *p)
{
  int i = 0;
  while ((fscanf(file, "%lf", p) != EOF) && (++i < n)) p++;
  return i;
}

void
print_usage(void)
{
  printf(
      "usage: a.out <omega> <zeta> - <beta> <gamma> [alpha]...\n"

  );
}

int
main(int argc, char **argv)
{
  // a.out period damp - beta gamma alpha
  double n = 6;
  double p[6];
  double dt = 0.1;

  read_load(stdin, n, p);

  int argi = 1;
  double omega  = atof(argv[argi++]);
  double zeta   = atof(argv[argi++]);
  argi++;
  conf.beta     = atof(argv[argi++]);
  conf.gamma    = atof(argv[argi++]);

  double M = 0.2533;
  double K = M*(omega*omega);
  // double C = 2.0*sqrt(K/M)*zeta;
  double C = 0.1592;
  double s = 1.0; // -386.4*M

  generalized_alpha(&conf,
    M, C, K, 
    s, n, p, dt);
}

