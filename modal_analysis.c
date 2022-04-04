#include <pthread.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <ctype.h>

#define handle_error_en(en, msg) \
       do { errno = en; perror(msg); exit(EXIT_FAILURE); } while (0)

struct load {
  int length;
  double dt;
  double series[];
};

void *thread_start(void *hist)
{
  /*
  for (n) {
    G[n] = L[n] / M[n];
    for (j)
      for (i)
        u[j][i] += shape[n][j]*D[n][i] * G[n];
  }
  */
  return 0;
}


struct mode {double *shape;};
struct step {double t, p, u, v, a;};

struct hist;

struct hist {
  int length, *past, *pres;
  const double *t, *p, *u, *v, *a;
  struct step *u_max, *v_max, *a_max;
  const void(*incr)(struct hist*);
  pthread_t thread_id;
};

void update(struct hist* hist)
{
  if (abs(hist->v[*hist->pres]) > hist->u_max->u) {
    hist->u_max->u = hist->u[*hist->pres];
    hist->u_max->v = hist->v[*hist->pres];
    hist->u_max->a = hist->a[*hist->pres];
  }
  if (abs(hist->v[*hist->pres]) > hist->v_max->u) {
    hist->v_max->u = hist->u[*hist->pres];
    hist->v_max->v = hist->v[*hist->pres];
    hist->v_max->a = hist->a[*hist->pres];
  }
  if (abs(hist->a[*hist->pres]) > hist->a_max->a) {
    hist->a_max->u = hist->u[*hist->pres];
    hist->a_max->v = hist->v[*hist->pres];
    hist->a_max->a = hist->a[*hist->pres];
  }
}

void rha_incr(struct hist* hist)
{
  update(hist);
  ++hist->past;
  ++hist->pres;
}

void rsa_incr(struct hist* hist)
{
  *hist->past = !(*hist->past);
  *hist->pres = !(*hist->pres);
}


void rha(int N, struct mode modes[N], struct load *load)
{
  struct step u_max[N], v_max[N], a_max[N];
  double u[N][2], v[N][2], a[N][2];
  struct hist history[N];
  for (int n = 0; n < N; n++) {
    history[n].length  = N,
    history[n].p       = load->series;
    history[n].u       = u[n];
    history[n].v       = u[n];
    history[n].a       = u[n];
    history[n].u_max   = &u_max[n];
    history[n].v_max   = &v_max[n];
    history[n].a_max   = &a_max[n];

    history[n].incr  = &rha_incr;
  }
}

void rsa(int N, struct mode modes[N], struct load *load)
{
  struct step u_max[N], v_max[N], a_max[N];

  struct hist history[N];
  for (int n = 0; n < N; n++) {
    history[n].length = load->length;
    history[n].p      = load->series;
    history[n].u      = NULL;
    history[n].v      = NULL;
    history[n].a      = NULL;
    history[n].u_max  = &u_max[n];
    history[n].v_max  = &v_max[n];
    history[n].a_max  = &a_max[n];

    history[n].incr  = &rsa_incr;
    int s = pthread_create(&history[n].thread_id, NULL, &thread_start, &history[n]);
    if (s != 0)
        handle_error_en(s, "pthread_create");
  }
}


