#include <tgmath.h>

void *newmark_step(
    double *u0, double *v0, double *p
    double dt,
    double gamma, double beta
    ) {
              khat = k + a1
              deltap = p[j] - p[j-1]
              deltaph = deltap + a2*v[j-1] + a3*a[j-1]
              deltau = deltaph/khat
              deltav = gamma*deltau/(beta*dtp) - gamma*v[j-1]/beta+dtp*(1-gamma/(2*beta))*a[j-1]
              deltaa = deltau/(beta*dtp**2) - v[j-1]/(beta*dtp) - a[j-1]/(2*beta)
              u[j] = u[j-1] + deltau
              v[j] = v[j-1] + deltav
              a[j] = a[j-1] + deltaa
}

void * _accel_spectrum(accel,dt,damping,per,gamma=1/2,beta=1/4,interp)
{
    pi2 = M_PI*M_PI;
    int numper = len(per);
    double m = 1.0;
    numdata = len(accel)
    t = np.arange(0,(numdata)*dt,dt);

    double u0 = 0.0, v0 = 0.0;
    SA = malloc((1+len(damping))*numper);
    SA[0,:] = per[:]

    for di,dmp in enumerate(damping) {
        for (int i = 0; i < numper; i++) {
            if (dt/per[i] > 0.02) {
                dtp = per[i]*0.02;
                dtpx = np.arange(0,max(t),dtp);
                dtpx = dtpx;
                accfrni = interp(t, accel)(dtpx);
                accfrn = accfrni[1:len(accfrni)-1];
                numdatan = len(accfrn)
                p  =  -m*accfrn;
            } else {
                dtp = dt;
                accfrn = accel;
                p = -m*accfrn;
                numdatan = numdata;
            }

            k = 4.0*pi2*m/(per[i]*per[i]);
            c = 2.0*dmp*sqrt(k/m);

            u = malloc((sizeof double)*numdatan);
            v = malloc((sizeof double)*numdatan);
            a = malloc((sizeof double)*numdatan);
            u[0] = u0
            v[0] = v0
            a[0] = (p[0] - c*v[0] - k*u[0])/m

            a1 = gamma*c/(beta*dtp) + m/(beta*dtp**2.0)
            a2 = m/(beta*dtp) + gamma*c/beta
            a3 = m/(2*beta) + dtp*(gamma/(2*beta)-1)*c
            khat = k + a1
            for (j = 1; j < numdatan; j++) {
                deltap = p[j] - p[j-1]
                deltaph = deltap + a2*v[j-1] + a3*a[j-1]
                deltau = deltaph/khat
                deltav = gamma*deltau/(beta*dtp) - gamma*v[j-1]/beta+dtp*(1-gamma/(2*beta))*a[j-1]
                deltaa = deltau/(beta*dtp**2) - v[j-1]/(beta*dtp) - a[j-1]/(2*beta)
                u[j] = u[j-1] + deltau
                v[j] = v[j-1] + deltav
                a[j] = a[j-1] + deltaa
            }
            atot = a + accfrn;
            SA[1+di,i] = abs(max(atot, key=abs))
        }
    }
    return SA;
}
