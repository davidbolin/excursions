

/* ------------------------ */
/* FROM GSL_MACHINE.H 		*/
/* ------------------------ */
/* Author:  B. Gough and G. Jungman */
#define GSL_DBL_EPSILON        2.2204460492503131e-16


/* ------------------------ */
/* FROM GSL_SYS.H 			*/
/* ------------------------ */

/* sys/gsl_sys.h
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

double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

/* ------------------------ */
/* FROM GSL_INLINE.H 		*/
/* ------------------------ */

/* gsl_inline.h
 * 
 * Copyright (C) 2008, 2009 Brian Gough
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
 
#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif
#define GSL_RANGE_COND(x) (x)

#endif 

/* ------------------------ */
/* FROM GSL_TYPES.H 		*/
/* ------------------------ */

/* gsl_types.h
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

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif


/* ------------------------ */
/* FROM GSL_CDF.H 			*/
/* ------------------------ */

/* cdf/gsl_cdf.h
 * 
 * Copyright (C) 2002 Jason H. Stover.
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

/* Author:  J. Stover */

#ifndef __GSL_CDF_H__
#define __GSL_CDF_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS          /* empty */
# define __END_DECLS            /* empty */
#endif

__BEGIN_DECLS 

double gsl_cdf_ugaussian_P (const double x);
double gsl_cdf_ugaussian_Q (const double x);

double gsl_cdf_ugaussian_Pinv (const double P);
double gsl_cdf_ugaussian_Qinv (const double Q);

double gsl_cdf_gaussian_P (const double x, const double sigma);
double gsl_cdf_gaussian_Q (const double x, const double sigma);

double gsl_cdf_gaussian_Pinv (const double P, const double sigma);
double gsl_cdf_gaussian_Qinv (const double Q, const double sigma);

__END_DECLS

#endif /* __GSL_CDF_H__ */
