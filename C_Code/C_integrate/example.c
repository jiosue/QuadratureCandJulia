#include <stdio.h>
#include "integrate.h"


double f(double x) {
    return x + 2 * x * x - 3 * x * x * x;
}


double df(double x) {
    // derivative of f
    return 1 + 4 * x - 9 * x * x;
}


int main() {
    double xa = -5; double xb = 10;
    int nfev = 100;

    printf("True integral: %.10f\n", f(xb) - f(xa));

    printf("Approximate with %d function evaluations\n", nfev);
    printf("\tRectangle rule: %.10f\n", Rectangular(df, xa, xb, nfev));
    printf("\tTrapezoid rule: %.10f\n", Trapezoidal(df, xa, xb, nfev));
    printf("\tSimpsons rule: %.10f\n\n", Simpson(df, xa, xb, nfev));

    return 0;
}
