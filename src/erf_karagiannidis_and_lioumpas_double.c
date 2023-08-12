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

/*  Coefficients from Karagiannidis and Lioumpas approximation.               */
#define A (-1.9800000000000000000000000000000000000000E+00)
#define B (+2.0117351207777607067259726964788258592758E+00)

/*  Approximation of erf(x) using the exponential function.                   */
double Erf_Double_Karagiannidis_and_Lioumpas(double x)
{
    /*  The function is odd in x so compute abs_x and use this.               */
    const double abs_x = tmpl_Double_Abs(x);

    /*  Compute the factors from K & L.                                       */
    const double exp_A_x = tmpl_Double_Exp(A * abs_x);
    const double exp_minus_x2 = tmpl_Double_Exp(-abs_x * abs_x);
    const double scale = exp_minus_x2 / (B * abs_x);
    const double erf_abs_x = 1.0 - (1.0 - exp_A_x) * scale;

    /*  The error function is odd. Return -erf_x for negative x.              */
    if (x < 0.0)
        return -erf_abs_x;

    /*  Otherwise we have a decent approximation for erf(x). Return this.     */
    return erf_abs_x;
}
/*  End of Erf_Double_Karagiannidis_and_Lioumpas.                             */

/*  Undefine everything in case someone wants to #include this file.          */
#undef A
#undef B
