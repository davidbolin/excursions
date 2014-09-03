#include <math.h>
#include "gsl_fix.h"
#include <stdlib.h>

/* ------------------------ */
/* FROM GSL_NAN.H			*/
/* ------------------------ */
/* gsl_nan.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif




/* ------------------------ */
/* FROM RATEVAL.H			*/
/* ------------------------ */

static double rat_eval (const double a[], const size_t na,
          const double b[], const size_t nb, const double x)
{
  size_t i, j;
  double u, v, r;

  u = a[na - 1];

  for (i = na - 1; i > 0; i--)
    {
      u = x * u + a[i - 1];
    }

  v = b[nb - 1];

  for (j = nb - 1; j > 0; j--)
    {
      v = x * v + b[j - 1];
    }

  r = u / v;

  return r;
}

/* ------------------------ */
/* FROM INFNAN.C			*/
/* ------------------------ */
/* sys/infnan.c
 * 
 * Copyright (C) 2001, 2004, 2007, 2010 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

double gsl_posinf (void)
{
  return gsl_fdiv (+1.0, 0.0);
}

double gsl_neginf (void)
{
  return gsl_fdiv (-1.0, 0.0);
}

/* ------------------------ */
/* FROM FDIV.C	 			*/
/* ------------------------ */
/* sys/fdiv.c
 * 
 * Copyright (C) 2001, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */


double gsl_fdiv (const double x, const double y)
{
  return x / y;
}

/* ------------------------ */
/* FROM GAUSS.C 			*/
/* ------------------------ */

/* cdf/gauss.c
 *
 * Copyright (C) 2002, 2004 Jason H. Stover.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * Computes the cumulative distribution function for the Gaussian
 * distribution using a rational function approximation.  The
 * computation is for the standard Normal distribution, i.e., mean 0
 * and standard deviation 1. If you want to compute Pr(X < t) for a
 * Gaussian random variable X with non-zero mean m and standard
 * deviation sd not equal to 1, find gsl_cdf_ugaussian ((t-m)/sd).
 * This approximation is accurate to at least double precision. The
 * accuracy was verified with a pari-gp script.  The largest error
 * found was about 1.4E-20. The coefficients were derived by Cody.
 *
 * References:
 *
 * W.J. Cody. "Rational Chebyshev Approximations for the Error
 * Function," Mathematics of Computation, v23 n107 1969, 631-637.
 *
 * W. Fraser, J.F Hart. "On the Computation of Rational Approximations
 * to Continuous Functions," Communications of the ACM, v5 1962.
 *
 * W.J. Kennedy Jr., J.E. Gentle. "Statistical Computing." Marcel Dekker. 1980.
 * 
 *  
 */
 
#ifndef M_1_SQRT2PI
#define M_1_SQRT2PI (M_2_SQRTPI * M_SQRT1_2 / 2.0)
#endif

#define SQRT32 (4.0 * M_SQRT2)
#define GAUSS_EPSILON  (GSL_DBL_EPSILON / 2)
#define GAUSS_XUPPER (8.572)
#define GAUSS_XLOWER (-37.519)
#define GAUSS_SCALE (16.0)

static double get_del (double x, double rational)
{
  double xsq = 0.0;
  double del = 0.0;
  double result = 0.0;

  xsq = floor (x * GAUSS_SCALE) / GAUSS_SCALE;
  del = (x - xsq) * (x + xsq);
  del *= 0.5;

  result = exp (-0.5 * xsq * xsq) * exp (-1.0 * del) * rational;

  return result;
}

/*
 * Normal cdf for fabs(x) < 0.66291
 */
static double
gauss_small (const double x)
{
  unsigned int i;
  double result = 0.0;
  double xsq;
  double xnum;
  double xden;

  const double a[5] = {
    2.2352520354606839287,
    161.02823106855587881,
    1067.6894854603709582,
    18154.981253343561249,
    0.065682337918207449113
  };
  const double b[4] = {
    47.20258190468824187,
    976.09855173777669322,
    10260.932208618978205,
    45507.789335026729956
  };

  xsq = x * x;
  xnum = a[4] * xsq;
  xden = xsq;

  for (i = 0; i < 3; i++)
    {
      xnum = (xnum + a[i]) * xsq;
      xden = (xden + b[i]) * xsq;
    }

  result = x * (xnum + a[3]) / (xden + b[3]);

  return result;
}

/*
 * Normal cdf for 0.66291 < fabs(x) < sqrt(32).
 */
static double
gauss_medium (const double x)
{
  unsigned int i;
  double temp = 0.0;
  double result = 0.0;
  double xnum;
  double xden;
  double absx;

  const double c[9] = {
    0.39894151208813466764,
    8.8831497943883759412,
    93.506656132177855979,
    597.27027639480026226,
    2494.5375852903726711,
    6848.1904505362823326,
    11602.651437647350124,
    9842.7148383839780218,
    1.0765576773720192317e-8
  };
  const double d[8] = {
    22.266688044328115691,
    235.38790178262499861,
    1519.377599407554805,
    6485.558298266760755,
    18615.571640885098091,
    34900.952721145977266,
    38912.003286093271411,
    19685.429676859990727
  };

  absx = fabs (x);

  xnum = c[8] * absx;
  xden = absx;

  for (i = 0; i < 7; i++)
    {
      xnum = (xnum + c[i]) * absx;
      xden = (xden + d[i]) * absx;
    }

  temp = (xnum + c[7]) / (xden + d[7]);

  result = get_del (x, temp);

  return result;
}

/*
 * Normal cdf for 
 * {sqrt(32) < x < GAUSS_XUPPER} union { GAUSS_XLOWER < x < -sqrt(32) }.
 */
static double
gauss_large (const double x)
{
  int i;
  double result;
  double xsq;
  double temp;
  double xnum;
  double xden;
  double absx;

  const double p[6] = {
    0.21589853405795699,
    0.1274011611602473639,
    0.022235277870649807,
    0.001421619193227893466,
    2.9112874951168792e-5,
    0.02307344176494017303
  };
  const double q[5] = {
    1.28426009614491121,
    0.468238212480865118,
    0.0659881378689285515,
    0.00378239633202758244,
    7.29751555083966205e-5
  };

  absx = fabs (x);
  xsq = 1.0 / (x * x);
  xnum = p[5] * xsq;
  xden = xsq;

  for (i = 0; i < 4; i++)
    {
      xnum = (xnum + p[i]) * xsq;
      xden = (xden + q[i]) * xsq;
    }

  temp = xsq * (xnum + p[4]) / (xden + q[4]);
  temp = (M_1_SQRT2PI - temp) / absx;

  result = get_del (x, temp);

  return result;
}

double
gsl_cdf_ugaussian_P (const double x)
{
  double result;
  double absx = fabs (x);

  if (absx < GAUSS_EPSILON)
    {
      result = 0.5;
      return result;
    }
  else if (absx < 0.66291)
    {
      result = 0.5 + gauss_small (x);
      return result;
    }
  else if (absx < SQRT32)
    {
      result = gauss_medium (x);

      if (x > 0.0)
        {
          result = 1.0 - result;
        }

      return result;
    }
  else if (x > GAUSS_XUPPER)
    {
      result = 1.0;
      return result;
    }
  else if (x < GAUSS_XLOWER)
    {
      result = 0.0;
      return result;
    }
  else
    {
      result = gauss_large (x);

      if (x > 0.0)
        {
          result = 1.0 - result;
        }
    }

  return result;
}

double
gsl_cdf_ugaussian_Q (const double x)
{
  double result;
  double absx = fabs (x);

  if (absx < GAUSS_EPSILON)
    {
      result = 0.5;
      return result;
    }
  else if (absx < 0.66291)
    {
      result = gauss_small (x);

      if (x < 0.0)
        {
          result = fabs (result) + 0.5;
        }
      else
        {
          result = 0.5 - result;
        }

      return result;
    }
  else if (absx < SQRT32)
    {
      result = gauss_medium (x);

      if (x < 0.0)
        {
          result = 1.0 - result;
        }

      return result;
    }
  else if (x > -(GAUSS_XLOWER))
    {
      result = 0.0;
      return result;
    }
  else if (x < -(GAUSS_XUPPER))
    {
      result = 1.0;
      return result;
    }
  else
    {
      result = gauss_large (x);

      if (x < 0.0)
        {
          result = 1.0 - result;
        }

    }

  return result;
}

double
gsl_cdf_gaussian_P (const double x, const double sigma)
{
  return gsl_cdf_ugaussian_P (x / sigma);
}

double
gsl_cdf_gaussian_Q (const double x, const double sigma)
{
  return gsl_cdf_ugaussian_Q (x / sigma);
}



/* ------------------------ */
/* FROM GAUSSINV.C 			*/
/* ------------------------ */
/* cdf/inverse_normal.c
 *
 * Copyright (C) 2002 Przemyslaw Sliwa and Jason H. Stover.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * Computes the inverse normal cumulative distribution function 
 * according to the algorithm shown in 
 *
 *      Wichura, M.J. (1988).
 *      Algorithm AS 241: The Percentage Points of the Normal Distribution.
 *      Applied Statistics, 37, 477-484.
 */



static double
small (double q)
{
  const double a[8] = { 3.387132872796366608, 133.14166789178437745,
    1971.5909503065514427, 13731.693765509461125,
    45921.953931549871457, 67265.770927008700853,
    33430.575583588128105, 2509.0809287301226727
  };

  const double b[8] = { 1.0, 42.313330701600911252,
    687.1870074920579083, 5394.1960214247511077,
    21213.794301586595867, 39307.89580009271061,
    28729.085735721942674, 5226.495278852854561
  };

  double r = 0.180625 - q * q;

  double x = q * rat_eval (a, 8, b, 8, r);

  return x;
}

static double
intermediate (double r)
{
  const double a[] = { 1.42343711074968357734, 4.6303378461565452959,
    5.7694972214606914055, 3.64784832476320460504,
    1.27045825245236838258, 0.24178072517745061177,
    0.0227238449892691845833, 7.7454501427834140764e-4
  };

  const double b[] = { 1.0, 2.05319162663775882187,
    1.6763848301838038494, 0.68976733498510000455,
    0.14810397642748007459, 0.0151986665636164571966,
    5.475938084995344946e-4, 1.05075007164441684324e-9
  };

  double x = rat_eval (a, 8, b, 8, (r - 1.6));

  return x;
}

static double
tail (double r)
{
  const double a[] = { 6.6579046435011037772, 5.4637849111641143699,
    1.7848265399172913358, 0.29656057182850489123,
    0.026532189526576123093, 0.0012426609473880784386,
    2.71155556874348757815e-5, 2.01033439929228813265e-7
  };

  const double b[] = { 1.0, 0.59983220655588793769,
    0.13692988092273580531, 0.0148753612908506148525,
    7.868691311456132591e-4, 1.8463183175100546818e-5,
    1.4215117583164458887e-7, 2.04426310338993978564e-15
  };

  double x = rat_eval (a, 8, b, 8, (r - 5.0));

  return x;
}

double
gsl_cdf_ugaussian_Pinv (const double P)
{
  double r, x, pp;

  double dP = P - 0.5;

  if (P == 1.0)
    {
      return GSL_POSINF;
    }
  else if (P == 0.0)
    {
      return GSL_NEGINF;
    }

  if (fabs (dP) <= 0.425)
    {
      x = small (dP);

      return x;
    }

  pp = (P < 0.5) ? P : 1.0 - P;

  r = sqrt (-log (pp));

  if (r <= 5.0)
    {
      x = intermediate (r);
    }
  else
    {
      x = tail (r);
    }

  if (P < 0.5)
    {
      return -x;
    }
  else
    {
      return x;
    }

}

double
gsl_cdf_ugaussian_Qinv (const double Q)
{
  double r, x, pp;

  double dQ = Q - 0.5;

  if (Q == 1.0)
    {
      return GSL_NEGINF;
    }
  else if (Q == 0.0)
    {
      return GSL_POSINF;
    }

  if (fabs (dQ) <= 0.425)
    {
      x = small (dQ);

      return -x;
    }

  pp = (Q < 0.5) ? Q : 1.0 - Q;

  r = sqrt (-log (pp));

  if (r <= 5.0)
    {
      x = intermediate (r);
    }
  else
    {
      x = tail (r);
    }

  if (Q < 0.5)
    {
      return x;
    }
  else
    {
      return -x;
    }
}


double
gsl_cdf_gaussian_Pinv (const double P, const double sigma)
{
  return sigma * gsl_cdf_ugaussian_Pinv (P);
}

double
gsl_cdf_gaussian_Qinv (const double Q, const double sigma)
{
  return sigma * gsl_cdf_ugaussian_Qinv (Q);
}
