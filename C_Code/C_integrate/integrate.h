#ifndef Integrate_H_INCLUDED
#define Integrate_H_INCLUDED

double Rectangular(double f(double), double xa, double xb, int nfev);

double Trapezoidal(double f(double), double xa, double xb, int nfev);

double Simpson(double f(double), double xa, double xb, int nfev);

#endif
