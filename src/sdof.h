// Claudio Perez
//
#ifndef SDOF_H
#define SDOF_H

struct sdof_alpha {
  double alpha_m,
         alpha_f,
         beta,
         gamma;
};

struct sdof_peaks {
    double max_displ,
           max_veloc,
           max_accel;
};

struct sdof_peaks
sdof_integrate_peaks_2(struct sdof_alpha* conf,
    double M, double C, double K,
    double scale, int n, const double *p, double dt);

#endif // SDOF_H
