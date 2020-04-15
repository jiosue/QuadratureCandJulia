#include "integrate.h"
#include <math.h>

// Numerical integration rules.


double Rectangular(double f(double), double xa, double xb, int nfev) {
    // Numerically integrate ``f`` from ``xa`` to ``xb`` with the
    // rectangular rule using ``nfev`` evaluations of ``f``.
    int M = nfev;
    double delta = (xb - xa) / (double)M;
    double sum = 0;
    for(int m=0; m<M; m++)
        sum += delta * f(xa + m * delta);
    return sum;
}


double Trapezoidal(double f(double), double xa, double xb, int nfev) {
    // Numerically integrate ``f`` from ``xa`` to ``xb`` with the
    // trapezoidal rule using ``nfev`` evaluations of ``f``.
    int M = nfev - 1;
    double delta = (xb - xa) / (double)M;
    double sum = (delta / 2) * (f(xa) + f(xb));
    sum += Rectangular(f, xa + delta, xb, M - 1);
    return sum;
}


double Simpson(double f(double), double xa, double xb, int nfev) {
    // Numerically integrate ``f`` from ``xa`` to ``xb`` with
    // Simpsons rule using ``nfev`` evaluations of ``f``.
    int M = (nfev - 1) / 2;
    double delta = (xb - xa) / (double)M;
    double sum = f(xa) + f(xb) + 4 * f(xa + delta / 2);
    for(int m=1; m<M; m++)
        sum += 2 * f(xa + (double)m * delta) +
               4 * f(xa + ((double)m + .5) * delta);
    return sum * delta / 6;
}
