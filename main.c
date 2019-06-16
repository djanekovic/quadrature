#include <stdio.h>
#include <math.h>

#include "quadrature.h"

double f_1(double x, double y, double z)
{
    return sin(x);
}

int main(int argc, char **argv)
{
    for (int i = 2; i <= 16; i++) {
        printf("%d: int sin(x)dx from 0 to pi: %.17g\n", i,
               compute_integral_line(f_1, 0, M_PI, i));
    }
}
