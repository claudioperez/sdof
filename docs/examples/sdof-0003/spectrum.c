#define WORK_SIZE 290

#include <sdof.h>
#include <stdio.h>
#include <stdlib.h>
extern struct sdof_alpha CONF;

static int
read_load(FILE* file, int n, double *p)
{
  int i = 0;
  while ((fscanf(file, "%lf", p) != EOF) && (++i < n)) p++;
  return i;
}

int main(int argc, char const *argv[]) {
  if (argc < 2) {
    printf("usage: spectrum <filename>\n");
    return -1;
  }

  FILE* f = fopen(argv[1], "r");

  if (f == 0) {
    printf("bad file: %s\n", argv[1]);
    return -1;
  }

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
  
  printf("# T \t Accel\n");
  for (int i=0; i< WORK_SIZE; i++) {
    double period = t_min + (t_max - t_min)/((double)WORK_SIZE)*((double)i);
    printf("%lf\t%lf\n", period, response[i].max_accel);
  }

  free(response);
}
