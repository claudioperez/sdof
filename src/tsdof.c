// Claudio Perez
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

#if defined(_WIN32)
#  include <Python.h>
   PyMODINIT_FUNC PyInit__tsdof(void) {}
#  define EXPORT __declspec(dllexport)

#elif defined(__EMSCRIPTEN__)
#  include <stdlib.h>
#  include <emscripten.h>
#  define EXPORT EMSCRIPTEN_KEEPALIVE

#else // *NIXs
#  define EXPORT
#endif


#ifndef M_PI
# define M_PI (3.14159265358979323846)
#endif
// Pre-define 2*pi
#define PI_PI (9.869604401089358)

#ifndef C11THREADS
# include <pthread.h>
#else
#include <threads.h>
# define pthread_create(a,b,c,d) thrd_create(a, c, d)
# define pthread_t               thrd_t
# define pthread_exit(a)         thrd_exit(a)
# define pthread_join(a,b)       thrd_join(a,b)
#endif

#define WORK_SIZE 290

// #ifndef NUM_THREADS
// # define NUM_THREADS 8
// #endif

struct sdof_peaks {
    double max_displ,
           max_veloc,
           max_accel;
};

extern struct sdof_alpha {
  double alpha_m,
         alpha_f,
         beta,
         gamma;
} CONF;


extern struct sdof_peaks
fsdof_peaks_2(struct sdof_alpha* conf,
    double M, double C, double K,
    double scale, int n, const double *p, double dt);

static int
read_load(FILE* file, int n, double *p)
{
  int i = 0;
  while ((fscanf(file, "%lf", p) != EOF) && (++i < n)) p++;
  return i;
}

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
//  printf("%d\t%lf\t%lf\t%lf\n", i+offset, period, C, K);
    resp[i*stride + offset] = fsdof_peaks_2(td.conf,   mass, C, K,
                                     scale, td.n, td.p, td.dt);
  }

  pthread_exit(EXIT_SUCCESS);
  return NULL;
}


EXPORT int
sdof_spectrum(struct sdof_alpha* conf,
              const double* load, const int n, const double dt, 
              const double t_min, const double t_max, const int n_periods,
              const double damp,
              int n_threads,
              struct sdof_peaks *response)
{
  pthread_t threads[n_threads];
  struct thread_data wkspace[n_threads];

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

  return 0;
}

#ifndef NO_MAIN
int main(int argc, char const *argv[]) {
  FILE* f = fopen("data/elCentro.txt", "r");

  double damp  = 0.02,
         t_min = 0.02,
         t_max = 3.00;

  double load[5000];
  int n = 1550;
  n = read_load(f, n, &load[0]);
  double dt = 0.02;

  const int n_threads = 8;

  struct sdof_peaks *response =(struct sdof_peaks *)
                               calloc(sizeof(struct sdof_peaks),WORK_SIZE);

  sdof_spectrum(&CONF, load, n, dt, t_min, t_max, WORK_SIZE, damp, n_threads, response);
  
  for (int i=0; i< WORK_SIZE; i++) {
    double period = t_min + (t_max - t_min)/((double)WORK_SIZE)*((double)i);
    printf("%lf\t%lf\n", period, response[i].max_accel);
  }

  free(response);
}
#endif

