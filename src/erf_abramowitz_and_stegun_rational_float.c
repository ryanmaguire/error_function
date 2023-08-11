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
#define A0 (1.0000000000F)
#define A1 (0.0705230784F)
#define A2 (0.0422820123F)
#define A3 (0.0092705272F)
#define A4 (0.0001520143F)
#define A5 (0.0002765672F)
#define A6 (0.0000430638F)

float Erf_Float_Abramowitz_and_Stegun_Rational(float x)
{
    float p, z;

    if (x > 6.0F)
        return 1.0F;
    else if (x < -6.0F)
        return - 1.0F;

    z = tmpl_Float_Abs(x);
    p = 1.0F / (A0 + z*(A1 + z*(A2 + z*(A3 + z*(A4 + z*(A5 + z*A6))))));
    p = p*p;
    p = p*p;
    p = p*p;
    p = p*p;

    if (x < 0.0F)
        return p - 1.0F;

    return 1.0F - p;
}
