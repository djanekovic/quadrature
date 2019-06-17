#include <stdio.h>
#include <math.h>

#include "quadrature.h"

double f_1(double x, double y, double z)
{
    return sin(x);
}

int main(int argc, char **argv)
{
    double exact = 2.0;
    quadrature_t q;

    for (int i = 2; i <= 16; i++) {
        init_quadrature(i, &q);
        printf("%d points: error = %.17g\n", i,
                exact - compute_integral_line(f_1, 0, M_PI, &q));
        destroy_quadrature(&q);
    }

    puts("");

    for (int i = 2; i <= 16; i++) {
        printf("%d points: error = %.17g\n", i,
                exact - compute_integral_line_alloc(f_1, 0, M_PI, i));
    }
}
