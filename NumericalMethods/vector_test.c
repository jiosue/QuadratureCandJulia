#include <stdio.h>
#include <stdlib.h>
#include "Vector.h"
/*
#define error_tol 0.0000001

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

Vector add(Vector a, Vector b) {
    Vector v;
    v.len = a.len;
    v.v = (double*)malloc(sizeof(double) * v.len);
    int i;
    for(i=0; i<v.len; i++) {
        v.v[i] = a.v[i] + b.v[i];
    }
    return v;    
}

Vector arb_add(int num, ...) {
    va_list valist; int i; 
    va_start(valist, num);
    Vector f = va_arg(valist, Vector);
    for(i=1; i<num; i++) {
        f = add(f, va_arg(valist, Vector));
    }
    va_end(valist);
    return f;
}

Vector subt(Vector a, Vector b) {
    //a - b;
    Vector v;
    v.len = a.len;
    v.v = (double*)malloc(sizeof(double) * v.len);
    int i;
    for(i=0; i<v.len; i++) {
        v.v[i] = a.v[i] - b.v[i];
    }
    return v; 
}

Vector scalar_mult(Vector a, double b) {
    Vector v;
    v.len = a.len;
    v.v = (double*)malloc(sizeof(double) * v.len);
    int i;
    for(i=0; i<v.len; i++) {
        v.v[i] = b*a.v[i];
    }
    return v;
}

double norm(Vector a) {
    double mag2 = 0; int i;
    for(i=0; i<a.len; i++) {
        mag2 += pow(a.v[i], 2);
    }
    return sqrt(mag2);
}

// Single Steps - Non-adaptive

Vector EulerStep(Vector f(double, Vector), double t, Vector y, double h) {
    return scalar_mult(f(t, y), h);
}

Vector ImprovedEulerStep(Vector f(double, Vector), double t, Vector y, double h) {
    Vector f0 = f(t, y);
    return scalar_mult(add(f0, f(t+h, add(y, scalar_mult(f0, h)))), h/2.0);
}

Vector MidpointStep(Vector f(double, Vector), double t, Vector y, double h) {
    return scalar_mult(f(t+h/2.0, add(y, scalar_mult(f(t, y), h/2.0))), h);
}

Vector RK4Step(Vector f(double, Vector), double t, Vector y, double h) {
    Vector f1 = f(t, y);
    Vector f2 = f(t+h/2.0, add(y, scalar_mult(f1, h/2.0)));
    Vector f3 = f(t+h/2.0, add(y, scalar_mult(f2, h/2.0)));
    Vector f4 = f(t+h, add(y, scalar_mult(f3, h)));
    return scalar_mult(arb_add(3, f1, f4, scalar_mult(add(f2, f3), 2)), h/6.0);
}

Data Stepper(Vector f(double, Vector), Vector method(Vector (double, Vector), double, Vector, double), double t0, Vector y0, double tf, double h) {
    Data d;
    d.len = (int)ceil((tf - t0) / h) + 1;
    d.ts = (double*)malloc(sizeof(double) * d.len);
    d.ys = (Vector*)malloc(sizeof(Vector) * d.len);
    
    int n;
    for(n=0; n<d.len; n++) {
        d.ts[n] = t0;
        d.ys[n] = y0;
        if(t0 + h > tf) {h = tf - t0;}
        y0 = add(y0, method(f, t0, y0, h));
        t0 += h;
    }
    return d;
}

// Single Steps - Adaptive

ErrorStep SimpleErrorStep(Vector f(double, Vector), double t, Vector y, double h) {
    ErrorStep e;
    e.step = RK4Step(f, t, y, h);
    e.error_est = norm(subt(e.step, MidpointStep(f, t, y, h)));
    return e;
}

ErrorStep DopriStep(Vector f(double, Vector), double t, Vector y, double h) {
    ErrorStep e; Vector f1, f2, f3, f4, f5, f6, f7;
    f1 = scalar_mult(f(t, y), h);
    f2 = scalar_mult(f(t + h/5.0, add(y, scalar_mult(f1, 1/5.0))), h);
    f3 = scalar_mult(f(t + 3*h/10.0, arb_add(3, y, scalar_mult(f1, 3/40.0), scalar_mult(f2, 9/40.0))), h);
    f4 = scalar_mult(f(t+4*h/5.0, arb_add(4, y, scalar_mult(f1, 44/45.0), scalar_mult(f2, -56/15.0), scalar_mult(f3, 32/9.0))), h);
    f5 = scalar_mult(f(t+8*h/9.0, arb_add(5, y, scalar_mult(f1, 19372/6561.0), scalar_mult(f2, -25360/2187.0), scalar_mult(f3, 64448/6561.0), scalar_mult(f4, -212/729.0))), h);
    f6 = scalar_mult(f(t+h, arb_add(6, y, scalar_mult(f1, 9017/3168.0), scalar_mult(f2, -355/33.0), scalar_mult(f3, 46732/5247.0), scalar_mult(f4, 49/176.0), scalar_mult(f5, -5103/18656.0))), h);
    e.step = arb_add(5, scalar_mult(f1, 35/384.0), scalar_mult(f3, 500/1113.0), scalar_mult(f4, 125/192.0), scalar_mult(f5, -2187/6784.0), scalar_mult(f6, 11/84.0));
    f7 = scalar_mult(f(t+h, add(y, e.step)), h);    
    e.error_est = norm(subt(e.step, arb_add(6, scalar_mult(f1, 5179/57600.0), scalar_mult(f3, 7571/16695.0), scalar_mult(f4, 393/640.0), scalar_mult(f5, -92097/339200.0), scalar_mult(f6, 187/2100.0), scalar_mult(f7, 1/40.0))));
    return e;
}


Data AdaptiveStepper(Vector f(double, Vector), ErrorStep method(Vector (double, Vector), double, Vector, double), double t0, Vector y0, double tf) {
    Data d;
    d.len = 1;
    d.ts = (double*)malloc(sizeof(double) * d.len);
    d.ys = (Vector*)malloc(sizeof(Vector) * d.len);
    d.ts[0] = t0; d.ys[0] = y0;
    
    double h = 0.001; double totalerror = 0;
    ErrorStep e;
    
    while(t0 < tf) {
        if(h < 0.00000001 && h < tf - t0) {
            printf("Step size effectively zero at t = %.10lf\n", t0);
            h = 0.001;
        } else {
            if(t0 + h > tf) {h = tf - t0;}
            e = method(f, t0, y0, h);
            if(e.error_est > error_tol) {h *= 0.75;}
            else {
                t0 += h; y0 = add(y0, e.step);

                //push t0 to ts and y0 to ys
                d.len++;
                d.ts = (double*)realloc(d.ts, sizeof(double) * d.len);
                d.ys = (Vector*)realloc(d.ys, sizeof(Vector) * d.len);
                d.ts[d.len-1] = t0; d.ys[d.len-1] = y0;
                
                if(e.error_est < error_tol / 10.0) {h *= 1.2;}
                totalerror += e.error_est;
            }
        }
    }
    printf("Estimated upper boud on total errorL %.10lf\n", totalerror);
    return d;
}

*/
/////////////////Test
Vector g(double t, Vector y) {
    //yddot = y;
    Vector v;
    v.len = 2;
    v.v = (double*)malloc(sizeof(double) * v.len);
    v.v[0] = y.v[1];
    v.v[1] = y.v[0];
    return v;
}
/*
void print_results(Vector method(Vector (double, Vector), double, Vector, double)) {
    double t0, tf, h;
    Vector y0;
    t0 = 0; tf = 1; h = 0.1;
    y0.len = 2;
    y0.v = (double*)malloc(sizeof(double) * y0.len);
    y0.v[0] = 1; y0.v[1] = 0;
    Data d = Stepper(g, method, t0, y0, tf, h);
    
    int i;
    for(i=0; i<d.len; i++) {
        printf("%f, %f\n", d.ts[i], d.ys[i].v[0]);
    }   
}

int main() {
    print_results(EulerStep);
    printf("\n\n");
    print_results(ImprovedEulerStep);
    printf("\n\n");
    print_results(MidpointStep);
    printf("\n\n");
    print_results(RK4Step);
    
	return 0;
}
*/

void print_results() {
    double t0, tf, h;
    Vector y0;
    t0 = 0; tf = 1; h = 0.1;
    y0.len = 2;
    y0.v = (double*)malloc(sizeof(double) * y0.len);
    y0.v[0] = 1; y0.v[1] = 0;
    Data d = AdaptiveStepper(g, DopriStep, t0, y0, tf);
    
    int i;
    for(i=0; i<d.len; i++) {
        printf("%f, %f\n", d.ts[i], d.ys[i].v[0]);
    }   
}

int main() {
    print_results();
    
	return 0;
}