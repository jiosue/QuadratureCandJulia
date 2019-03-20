#ifndef Vector_H_INCLUDED
#define Vector_H_INCLUDED

typedef struct Vector {
    int len;
    double *v;
} Vector;

typedef struct Data {
    int len;
    double *ts;
    Vector *ys;
} Data;

typedef struct ErrorStep {
    Vector step;
    double error_est;
} ErrorStep;

//  Vector Methods

Vector add(Vector a, Vector b);

Vector arb_add(int num, ...);

Vector subt(Vector a, Vector b);

Vector scalar_mult(Vector a, double b);

double norm(Vector a);

// Single Steps - Non-adaptive

Vector EulerStep(Vector f(double, Vector), double t, Vector y, double h);

Vector ImprovedEulerStep(Vector f(double, Vector), double t, Vector y, double h);

Vector MidpointStep(Vector f(double, Vector), double t, Vector y, double h);

Vector RK4Step(Vector f(double, Vector), double t, Vector y, double h);

Data Stepper(Vector f(double, Vector), Vector method(Vector (double, Vector), double, Vector, double), double t0, Vector y0, double tf, double h);

// Single Steps - Adaptive

ErrorStep SimpleErrorStep(Vector f(double, Vector), double t, Vector y, double h);

ErrorStep DopriStep(Vector f(double, Vector), double t, Vector y, double h);

Data AdaptiveStepper(Vector f(double, Vector), ErrorStep method(Vector (double, Vector), double, Vector, double), double t0, Vector y0, double tf);

#endif // Vector_H_INCLUDED
