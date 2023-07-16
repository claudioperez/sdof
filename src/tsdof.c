#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
# define M_PI (3.14159265358979323846)
#endif
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

#define WORK_SIZE 200

#ifndef NUM_THREADS
# define NUM_THREADS 8
#endif

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


struct sdof_peaks
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

  const double scale = 1.0,
               mass  = 1.0,
               damp  = td.damp;

  for (int offset = 0 ; offset < stride; offset++) {
    double period = td.t_min + td.t_slp*((double)(i+offset));
    double K = 4.0*PI_PI*mass/(period*period);
    double C = 4.0*damp*M_PI/period;
//  printf("%d\t%lf\t%lf\t%lf\n", i+offset, period, C, K);
    resp[i + offset] = fsdof_peaks_2(td.conf,   mass, C, K,
                                     scale, td.n, td.p, td.dt);
  }

//thrd_exit(EXIT_SUCCESS);
  return NULL;
}


int
sdof_spectrum(const double* load, const int n, const double dt, 
              const double t_min, const double t_max, const int n_periods,
              const double damp,
              struct sdof_peaks *response)
{
  pthread_t threads[NUM_THREADS];
  struct thread_data wkspace[NUM_THREADS];

  double slope = (t_max - t_min)/((double)n_periods);

  for (int i = 0; i < NUM_THREADS; i++) {
    wkspace[i].damp   =  damp; //0.1592;
    wkspace[i].n      =     n;
    wkspace[i].dt     =    dt;
    wkspace[i].p      =  load;
    wkspace[i].conf   = &CONF;
    wkspace[i].stride = n_periods/NUM_THREADS;

    wkspace[i].t_min  = t_min; // min. period
    wkspace[i].t_slp  = slope;
    wkspace[i].response = response;
    wkspace[i].thread_index = i*(n_periods/NUM_THREADS);

    pthread_create(&threads[i], NULL, &run_peaks, (void *)&wkspace[i]);
  }

  for(int i = 0; i < NUM_THREADS; i++)
    pthread_join(threads[i], NULL);

  return 0;
}

int main(int argc, char const *argv[]) {
  FILE* f = fopen("data/elCentro.txt", "r");

  double damp  = 0.02,
         t_min = 0.02,
         t_max = 3.00;

  double load[5000];
  int n = 1550;
  n = read_load(f, n, &load[0]);
  double dt = 0.02;

  struct sdof_peaks *response =(struct sdof_peaks *)
                               calloc(sizeof(struct sdof_peaks),WORK_SIZE);

  sdof_spectrum(load, n, dt, t_min, t_max, WORK_SIZE, damp, response);
  
  for (int i=0; i< WORK_SIZE; i++) {
    double period = t_min + (t_max - t_min)/((double)WORK_SIZE)*((double)i);
    printf("%lf\t%lf\n", period, response[i].max_accel);
  }

  free(response);
}

