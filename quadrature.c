#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "quadrature.h"

/* Compute line integral using n-points quadrature */
double compute_integral_line(double (*f)(double x, double y, double z),
                             double a, double b, int n)
{
    double *tmp, *pw, sum = 0;
    int size, num_points;

    double A = (b - a)/2.0;
    double B = (a + b)/2.0;

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
            fprintf(stderr, "Highest supported order is 8, defaulting to 8.\n");
            n = 16;
        case 16:
            size = QUAD_1D_16_LEN;
            tmp = (double[QUAD_1D_16_LEN]){ QUAD_1D_16 };
            break;
    }
    pw = malloc(size*sizeof(double));
    memcpy(pw, tmp, size * sizeof(double));

    num_points = size/2;


    if (n % 2 == 1) {
        sum = pw[1] * f(B, 0, 0);

        for (int i = 1; i < num_points; i++) {
            int _i2 = i * 2;
            double _f = f(A * pw[_i2] + B, 0, 0) + f(B - pw[_i2] * A, 0, 0);
            sum += pw[_i2+1] * _f;
        }
    } else {
        for (int i = 0; i < num_points; i++) {
            int _i2 = i * 2;
            double _f = f(A * pw[_i2] + B, 0, 0) + f(B - pw[_i2] * A, 0, 0);
            sum += pw[_i2+1] * _f;
        }
    }

    return A*sum;
}
