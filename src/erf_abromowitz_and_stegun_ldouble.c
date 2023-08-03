/******************************************************************************
 *                                  LICENSE                                   *
 ******************************************************************************
 *  This file is part of error_function.                                      *
 *                                                                            *
 *  error_function is free software: you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  error_function is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with error_function.  If not, see <https://www.gnu.org/licenses/>.  *
 ******************************************************************************/

/*  Absolute value and exponential function provided here.                    */
#include <libtmpl/include/tmpl_math.h>

/*  Function prototype provided here.                                         */
#include "error_function.h"

/*  This method of evaluation comes straight from Abramowitz and Stegun's     *
 *  textbook. They write that the maximum absolute error is bounded by 10^-7  *
 *  which makes this method appropriate for single precision, but nothing     *
 *  higher (double, extended, etc.).                                          */
#define A0 (+2.548295920E-01L)
#define A1 (-2.844967360E-01L)
#define A2 (+1.421413741E+00L)
#define A3 (-1.453152027E+00L)
#define A4 (+1.061405429E+00L)
#define P0 (+3.275911000E-01L)

/*  Polynomial evaluation via Horner's method.                                */
#define POLY_EVAL(z) z*(A0 + z*(A1 + z*(A2 + z*(A3 + z*A4))))

/*  The expansion is performed in terms of a rational function of z.          */
#define RATIONAL_TRANSFORM(z) 1.0L / (1.0L + P0*z)

/*  Rational approximation of erf(x) using the exponential function.          */
long double Erf_LDouble_Abramowitz_and_Stegun(long double x)
{
    /*  The function is odd in x so compute abs_x and use this.               */
    const long double abs_x = tmpl_LDouble_Abs(x);

    /*  Following Abramowitz and Stegun, we first compute exp(-x^2).          */
    const long double exp_minus_x_sq = tmpl_LDouble_Exp(-x*x);

    /*  Transform the variable to 1 / (1 + px), p is given above.             */
    const long double z = RATIONAL_TRANSFORM(abs_x);

    /*  Evaluate the polynomial using Horner's method.                        */
    const long double poly = POLY_EVAL(z);

    /*  The error function can be approximated as follows.                    */
    const long double erf_x = 1.0L - poly*exp_minus_x_sq;

    /*  The error function is odd. Return -erf_x for negative x.              */
    if (x < 0.0L)
        return -erf_x;

    /*  Otherwise we have a decent approximation for erf(x). Return this.     */
    return erf_x;
}
/*  End of Erf_LDouble_Abramowitz_and_Stegun.                                 */

/*  Undefine everything in case someone wants to #include this file.          */
#undef A0
#undef A1
#undef A2
#undef A3
#undef A4
#undef POLY_EVAL
#undef RATIONAL_TRANSFORM
