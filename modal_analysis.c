"modal_history "


/*

    w1   w2
   o11  o21
   o12  o22
   o13  o23
 */
#define handle_error_en(en, msg) \
       do { errno = en; perror(msg); exit(EXIT_FAILURE); } while (0)

void modal_compose(int i, struct hist[i])
{
  for (n)
    G[n] = L[n] / M[n];
    for (j)
      for (i)
        u[j][i] += shape[n][j]*D[n][i] * G[n]
}

struct step {double t, p, u, v, a;};

struct hist;
struct hist {
  int size, *past, *pres;
  const double *t, *p, *u, *v, *a;
  struct step *u_max, *v_max, *a_max;
  const void(*incr)(struct hist*);
  pthread_t thread_id;
};

void rsa_incr(struct hist* hist)
{
  hist->past = !hist->past;
  hist->pres = !hist->pres;
}

void update(struct hist* hist)
{
  if (abs(hist->v[*hist->pres]) > *hist->u_max->u) {
    *hist->u_max->u = hist->u[hist->pres];
    *hist->u_max->v = hist->v[hist->pres];
    *hist->u_max->a = hist->a[hist->pres];
  }
  if (abs(hist->v[*hist->pres]) > *hist->v_max->u) {
    *hist->v_max->u = hist->u[hist->pres];
    *hist->v_max->v = hist->v[hist->pres];
    *hist->v_max->a = hist->a[hist->pres];
  }
  if (abs(hist->a[*hist->pres]) > *hist->a_max->a) {
    *hist->a_max->u = hist->u[hist->pres];
    *hist->a_max->v = hist->v[hist->pres];
    *hist->a_max->a = hist->a[hist->pres];
  }
}

void rha_incr(struct hist* hist)
{
  rha_update(hist);
  ++hist->past;
  ++hist->pres;
}

void rha(int N, struct mode modes[N], struct load *load)
{
  struct step u_max[N], v_max[N], a_max[N];
  double u[N][2], v[N][2], a[N][2];
  for (n = 0; n < N; n++) {
    history[n].size  = N,
    history[n].p     = load->series;
    history[n].u     = u[n];
    history[n].v     = u[n];
    history[n].a     = u[n];
    history[n].u_max = &u_max[n];
    history[n].v_max = &v_max[n];
    history[n].a_max = &a_max[n];

    history[n].incr  = &rha_incr;
    pthread_create(&thread_id, );
  }
}

void rsa(int N, struct mode modes[N], struct load *load)
{
  struct step u_max[N], v_max[N], a_max[N];

  struct hist history[N];
  for (n = 0; n < N; n++) {
    history[n].length = load->length;
    history[n].p      = load->series;
    history[n].u      = NULL;
    history[n].v      = NULL;
    history[n].a      = NULL;
    history[n].u_max  = &u_max[n];
    history[n].v_max  = &v_max[n];
    history[n].a_max  = &a_max[n];

    history[n].incr  = &rsa_incr;
    s = pthread_create(&history[n].thread_id, NULL, &thread_start, &history[n]);
    if (s != 0)
        handle_error_en(s, "pthread_create");
  }
}


