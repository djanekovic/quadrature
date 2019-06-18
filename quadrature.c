#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "quadrature.h"

int _init_quadrature_1(int n, quadrature_t *q)
{
    int size;
    double *tmp;

    switch(n) {
        case 2:
            size = QUAD_1D_2_LEN;
            tmp = (double[QUAD_1D_2_LEN]){ QUAD_1D_2 };
            break;
        case 3:
            size = QUAD_1D_3_LEN;
            tmp = (double[QUAD_1D_3_LEN]){ QUAD_1D_3 };
            break;
        case 4:
            size = QUAD_1D_4_LEN;
            tmp = (double[QUAD_1D_4_LEN]){ QUAD_1D_4 };
            break;
        case 5:
            size = QUAD_1D_5_LEN;
            tmp = (double[QUAD_1D_5_LEN]){ QUAD_1D_5 };
            break;
        case 6:
            size = QUAD_1D_6_LEN;
            tmp = (double[QUAD_1D_6_LEN]){ QUAD_1D_6 };
            break;
        case 7:
            size = QUAD_1D_7_LEN;
            tmp = (double[QUAD_1D_7_LEN]){ QUAD_1D_7 };
            break;
        case 8:
            size = QUAD_1D_8_LEN;
            tmp = (double[QUAD_1D_8_LEN]){ QUAD_1D_8 };
            break;
        case 9:
            size = QUAD_1D_9_LEN;
            tmp = (double[QUAD_1D_9_LEN]){ QUAD_1D_9 };
            break;
        case 10:
            size = QUAD_1D_10_LEN;
            tmp = (double[QUAD_1D_10_LEN]){ QUAD_1D_10 };
            break;
        case 11:
            size = QUAD_1D_11_LEN;
            tmp = (double[QUAD_1D_11_LEN]){ QUAD_1D_11 };
            break;
        case 12:
            size = QUAD_1D_12_LEN;
            tmp = (double[QUAD_1D_12_LEN]){ QUAD_1D_12 };
            break;
        case 13:
            size = QUAD_1D_13_LEN;
            tmp = (double[QUAD_1D_13_LEN]){ QUAD_1D_13 };
            break;
        case 14:
            size = QUAD_1D_14_LEN;
            tmp = (double[QUAD_1D_14_LEN]){ QUAD_1D_14 };
            break;
        case 15:
            size = QUAD_1D_15_LEN;
            tmp = (double[QUAD_1D_15_LEN]){ QUAD_1D_15 };
            break;
        default:
            fprintf(stderr, "Order %d unsupported, defaulting to 16.\n", n);
            n = 16;
        case 16:
            size = QUAD_1D_16_LEN;
            tmp = (double[QUAD_1D_16_LEN]){ QUAD_1D_16 };
            break;
    }

    q->pw = malloc(size * sizeof(double));
    if (q->pw == NULL) {
        fprintf(stderr, "Allocating memory for point-weight failed\n");
        return 1;
    }

    memcpy(q->pw, tmp, size * sizeof(double));
    q->size = size;
    q->n = n;

    return 0;
}

int _init_quadrature_2(int n, quadrature_t *q)
{
    double *tmp;
    int size;

    switch (n) {
        case 1:
            tmp = (double[QUAD_2D_1_LEN]){ QUAD_2D_1 };
            size = QUAD_2D_1_LEN;
            break;
        case 2:
            tmp = (double[QUAD_2D_2_LEN]){ QUAD_2D_2 };
            size = QUAD_2D_2_LEN;
            break;
        case 3:
            tmp = (double[QUAD_2D_3_LEN]){ QUAD_2D_3 };
            size = QUAD_2D_3_LEN;
            break;
        case 4:
            tmp = (double[QUAD_2D_4_LEN]){ QUAD_2D_4 };
            size = QUAD_2D_4_LEN;
            break;
        case 5:
            tmp = (double[QUAD_2D_5_LEN]){ QUAD_2D_5 };
            size = QUAD_2D_5_LEN;
            break;
        case 6:
            tmp = (double[QUAD_2D_6_LEN]){ QUAD_2D_6 };
            size = QUAD_2D_6_LEN;
            break;
        case 7:
            tmp = (double[QUAD_2D_7_LEN]){ QUAD_2D_7 };
            size = QUAD_2D_7_LEN;
            break;
        default:
            fprintf(stderr, "Order %d unsupported, defaulting to 8.\n", n);
            n = 8;
        case 8:
            tmp = (double[QUAD_2D_8_LEN]){ QUAD_2D_8 };
            size = QUAD_2D_8_LEN;
            break;
    }

    q->pw = malloc(size * sizeof(double));
    if (q->pw == NULL) {
        fprintf(stderr, "Allocating memory for point-weight failed\n");
        return 1;
    }

    memcpy(q->pw, tmp, size * sizeof(double));
    q->size = size/3;
    q->n = n;

    return 0;
}

int init_quadrature(int dimensions, int n, quadrature_t *q)
{
    int retval;

    switch(dimensions) {
        case 1:
            retval = _init_quadrature_1(n, q);
            break;
        case 2:
            retval = _init_quadrature_2(n, q);
            break;
        default:
            fprintf(stderr, "Dimension %d unsupported\n", dimensions);
            retval = 1;
    }
    return retval;
}

void destroy_quadrature(quadrature_t *q)
{
    free(q->pw);
}

double compute_integral_line(double (*f)(double x, double y, double z),
                             double a, double b, quadrature_t *q)
{
    double A, B, sum = 0;
    int num_points;

    A = (b - a)/2.0;
    B = (a + b)/2.0;

    num_points = q->size>>1;

    /* if number of points is odd */
    if (q->n & 1) {
        sum = q->pw[1] * f(B, 0, 0);

        for (int i = 1; i < num_points; i++) {
            int _i2 = i << 1;
            double _f = f(A * q->pw[_i2] + B, 0, 0) + f(B - q->pw[_i2] * A, 0, 0);
            sum += q->pw[_i2+1] * _f;
        }
    } else {
        for (int i = 0; i < num_points; i++) {
            int _i2 = i << 1;
            double _f = f(A * q->pw[_i2] + B, 0, 0) + f(B - q->pw[_i2] * A, 0, 0);
            sum += q->pw[_i2+1] * _f;
        }
    }

    return A*sum;
}

/* Compute line integral using n-points quadrature */
double compute_integral_line_alloc(double (*f)(double x, double y, double z),
                                   double a, double b, int n)
{
    double A, B, sum = 0;
    int num_points;
    quadrature_t q;

    A = (b - a)/2.0;
    B = (a + b)/2.0;

    if(_init_quadrature_1(n, &q)) {
        fprintf(stderr, "Quadrature table init failed, exiting...\n");
        exit(1);
    }

    num_points = q.size>>1;

    /* if number of points is odd */
    if (n & 1) {
        sum = q.pw[1] * f(B, 0, 0);

        for (int i = 1; i < num_points; i++) {
            int _i2 = i << 1;
            double _f = f(A * q.pw[_i2] + B, 0, 0) + f(B - q.pw[_i2] * A, 0, 0);
            sum += q.pw[_i2+1] * _f;
        }
    } else {
        for (int i = 0; i < num_points; i++) {
            int _i2 = i << 1;
            double _f = f(A * q.pw[_i2] + B, 0, 0) + f(B - q.pw[_i2] * A, 0, 0);
            sum += q.pw[_i2+1] * _f;
        }
    }

    destroy_quadrature(&q);

    return A*sum;
}

double compute_integral_triangle_alloc(double (*f)(double x, double y, double z),
                                       int n)
{
    double sum = 0;
    quadrature_t q;

    if (_init_quadrature_2(n, &q)) {
        fprintf(stderr, "Quadrature table init failed, exiting...\n");
        exit(1);
    }

    for (int i = 0; i < q.size; i++) {
        int _i3 = i * 3;

        sum += q.pw[_i3+2]*f(q.pw[_i3], q.pw[_i3+1], 0);
    }

    destroy_quadrature(&q);

    return 0.5 * sum;
}

double compute_integral_triangle(double (*f)(double x, double y, double z),
                                 quadrature_t *q)
{
    double sum = 0;

    for (int i = 0; i < q->size; i++) {
        int _i3 = i * 3;

        sum += q->pw[_i3+2]*f(q->pw[_i3], q->pw[_i3+1], 0);
    }

    return 0.5 * sum;
}
