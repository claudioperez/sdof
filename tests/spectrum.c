#define WORK_SIZE 290
#include <stdio.h>
extern struct sdof_alpha CONF;

static int
read_load(FILE* file, int n, double *p)
{
  int i = 0;
  while ((fscanf(file, "%lf", p) != EOF) && (++i < n)) p++;
  return i;
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
