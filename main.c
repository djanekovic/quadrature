#include <stdio.h>
#include <math.h>

#include "quadrature.h"

double f_1(double x, double y, double z)
{
    return 2.0;
}

int main(int argc, char **argv)
{
    double exact = 2.0;
    quadrature_t q;

    for (int i = 2; i <= 16; i++) {
        init_quadrature(1, i, &q);
        printf("%d points: error = %.17g\n", i,
                exact - compute_integral_line(f_1, 0, 2.0, &q));
        destroy_quadrature(&q);
    }

    for (int i = 2; i <= 16; i++) {
        printf("%d points: error = %.17g\n", i,
                M_PI * 2.0 - compute_integral_line_alloc(f_1, 0, M_PI, i));
    }

    for (int i = 1; i <= 8; i++) {
        init_quadrature(2, i, &q);
        printf("%d points: error = %.17g\n", i,
                0.5*2 - compute_integral_triangle_alloc(f_1, i));
        destroy_quadrature(&q);
    }
}
