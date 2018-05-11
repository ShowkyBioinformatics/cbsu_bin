/* tajD.c -- Calculates Tajima's D (1989) and Fu and Li's D* (1993).
   Freq. spectrum etc. are embedded in the source code.  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
  int a, b, N, S, sing;
  double a1, a2, b1, b2, c1, c2, e1, e2, c, d, mu, nu, k, num, den, D, sqrt();
  double num2, den2, Dstar, pi, theta;

  if (argc==1) {
    printf("Usage: sample size segsites pi_total singletons\n"); exit(0);}

  N = atoi(argv[1]);
  S = atoi(argv[2]);
  k = atof(argv[3]);
  sing = atoi(argv[4]);

  a1 = a2 = 0.0;
  for (b=1; b<N; ++b) {
    a1 += 1. / (double) b;
    a2 += 1. / (double) (b * b);
  }
  b1 = (double) (N + 1.) / (3. * (N - 1.));
  b2 = (double) (2. * N * N + 2. * N + 6.) / (9. * N * (N - 1.));
  c1 = b1 - 1. / a1;
  c2 = b2 - (double) (N + 2.) / (a1 * N) + a2 / (a1 * a1);
  e1 = c1 / a1;
  e2 = c2 / (a1 * a1 + a2);
  num = k - (double) S / a1;
  den = sqrt (e1 * S + e2 * S * (S - 1));
  D = num / den;

  c = (double) (2. * (N * a1 - 2. * (N - 1.))) / ((N - 1.) * (N - 2.));
  d = (double) (1.5 - (2. * (a1 + 1. / N) - 3.) / (N - 2.) - 1. / N);
  d *= 2. / (N - 1.);
  d += c + (N - 2.) / (N * N - 2. * N + 1.);
  nu = (a2 * N * N) / ((N - 1.) * (N - 1.)) + a1 * a1 * d;
  nu -= (2. * N * a1 * (a1 +1)) / ((N - 1.) * (N - 1.));
  nu /= (a1 * a1 + a2);
  mu = N * (a1 - N / (N - 1.)) / (N - 1.) - nu;

  num2 = (double) N * S / (N - 1.) - a1 * sing;
  den2 = sqrt (mu * S + nu * S * S);
  Dstar = num2 / den2;


  theta = (double) S / a1;

  printf("Tajima's D = %f\nFu and Li's D* = %f\n\n", D, Dstar);
}  
