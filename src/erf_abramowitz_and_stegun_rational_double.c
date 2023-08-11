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

/*  Absolute value function provided here.                                    */
#include <libtmpl/include/tmpl_math.h>

/*  Function prototype provided here.                                         */
#include "error_function.h"

/*  Coefficients from Abramowitz and Stegun.                                  */
#define A0 (1.0000000000)
#define A1 (0.0705230784)
#define A2 (0.0422820123)
#define A3 (0.0092705272)
#define A4 (0.0001520143)
#define A5 (0.0002765672)
#define A6 (0.0000430638)

double Erf_Double_Abramowitz_and_Stegun_Rational(double x)
{
    double p, z;

    if (x > 6.0)
        return 1.0;
    else if (x < -6.0)
        return - 1.0;

    z = tmpl_Double_Abs(x);
    p = 1.0 / (A0 + z*(A1 + z*(A2 + z*(A3 + z*(A4 + z*(A5 + z*A6))))));
    p = p*p;
    p = p*p;
    p = p*p;
    p = p*p;

    if (x < 0.0)
        return p - 1.0;

    return 1.0 - p;
}
