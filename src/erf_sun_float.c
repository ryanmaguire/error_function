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

/*  Square root and exponential functions provided here.                      */
#include <libtmpl/include/tmpl_math.h>

/*  Function prototype provided here.                                         */
#include "error_function.h"

/*  Shift coefficients for the rational approximation.                        */
#define C0 (+8.42697144e-01F)

/*  Scale factor for the leading term of the Taylor series.                   */
#define TWO_BY_ROOT_PI (+1.1283791670955125738961589031215451716881012586E+00F)

/*  Coefficients for the rational approximation to erf on [0, 0.84375].       */
#define PA0 (+1.28379166E-01F)
#define PA1 (-3.36030394E-01F)
#define PA2 (-1.86260219E-03F)

#define QA0 (+1.00000000E+00F)
#define QA1 (+3.12324286E-01F)
#define QA2 (+2.16070302E-02F)
#define QA3 (-1.98859419E-03F)

/*  Helper macros for evaluating polynomials using Horner's method.           */
#define POLYPA_EVAL(t) PA0 + t*(PA1 + t*PA2)
#define POLYQA_EVAL(t) QA0 + t*(QA1 + t*(QA2 + t*QA3))

/*  Coefficients for the rational approximation to erf on [0.84375, 1.25].    */
#define PB0 (+3.64939137E-06F)
#define PB1 (+4.15109694E-01F)
#define PB2 (-1.65179938E-01F)
#define PB3 (+1.10914491E-01F)

#define QB0 (+1.00000000E+00F)
#define QB1 (+6.02074385E-01F)
#define QB2 (+5.35934687E-01F)
#define QB3 (+1.68576106E-01F)
#define QB4 (+5.62181212E-02F)

/*  Helper macros for evaluating polynomials using Horner's method.           */
#define POLYPB_EVAL(t) PB0 + t*(PB1 + t*(PB2 + t*PB3))
#define POLYQB_EVAL(t) QB0 + t*(QB1 + t*(QB2 + t*(QB3 + t*QB4)))

/*  Coefficients for the exponential approximation to erf on [1.25, 2.8571].  */
#define PC0 (-9.87132732E-03F)
#define PC1 (-5.53605914E-01F)
#define PC2 (-2.17589188E+00F)
#define PC3 (-1.43268085E+00F)

#define QC0 (+1.00000000E+00F)
#define QC1 (+5.45995426E+00F)
#define QC2 (+6.69798088E+00F)
#define QC3 (+1.43113089E+00F)
#define QC4 (-5.77397496E-02F)

/*  Helper macros for evaluating polynomials using Horner's method.           */
#define POLYPC_EVAL(t) PC0 + t*(PC1 + t*(PC2 + t*PC3))
#define POLYQC_EVAL(t) QC0 + t*(QC1 + t*(QC2 + t*(QC3 + t*QC4)))

/*  Coefficients for the exponential approximation to erf on [2.8571, 4.0].   */
#define PD0 (-9.86494310E-03F)
#define PD1 (-6.25171244E-01F)
#define PD2 (-6.16498327E+00F)
#define PD3 (-1.66696873E+01F)
#define PD4 (-9.53764343E+00F)

#define QD0 (+1.00000000E+00F)
#define QD1 (+1.26884899E+01F)
#define QD2 (+4.51839523E+01F)
#define QD3 (+4.72810211E+01F)
#define QD4 (+8.93033314E+00F)

/*  Helper macros for evaluating polynomials using Horner's method.           */
#define POLYPD_EVAL(t) PD0 + t*(PD1 + t*(PD2 + t*(PD3 + t*PD4)))
#define POLYQD_EVAL(t) QD0 + t*(QD1 + t*(QD2 + t*(QD3 + t*QD4)))

/*  Computes erf(x) using the algorithm outlined in the comments of SunOS's   *
 *  implementation of libm. This does not use their actual code, which makes  *
 *  use of macros assuming IEEE-754 support. Instead this uses the same idea  *
 *  but with code that is compliant with the C89 standard and quite portable. */
float Erf_Float_Sun(float x)
{
    /*  erf(x) is odd. We'll compute erf(|x|) and check the sign of x later.  */
    float erf_abs_x;
    const float abs_x = tmpl_Float_Abs(x);

    /*  For |x| >= 4 we have erf(|x|) = 1 to single precision. Return this.   */
    if (abs_x >= 4.0F)
        erf_abs_x = 1.0F;

    /*  For very small values, use the leading term of the Taylor series.     */
    else if (abs_x < 6.103515625e-05F)
        return TWO_BY_ROOT_PI*x;

    /*  Slightly larger, use a Remez rational minimax approximation.          */
    else if (abs_x < 0.84375F)
    {
        const float z = abs_x * abs_x;
        const float P = POLYPA_EVAL(z);
        const float Q = POLYQA_EVAL(z);
        return x + x * P/Q;
    }

    /*  For larger values, use a rational approximation centered at x = 1.    */
    else if (abs_x < 1.25F)
    {
        const float z = abs_x - 1.0F;
        const float P = POLYPB_EVAL(z);
        const float Q = POLYQB_EVAL(z);
        erf_abs_x = C0 + P/Q;
    }

    /*  And for all other values, use an exponentiated rational approximation.*/
    else
    {
        const float rcpr_abs_x = 1.0F / abs_x;
        const float z = rcpr_abs_x * rcpr_abs_x;
        float P, Q, factor;

        if (abs_x < 2.857142857142857F)
        {
            P = POLYPC_EVAL(z);
            Q = POLYQC_EVAL(z);
        }

        else
        {
            P = POLYPD_EVAL(z);
            Q = POLYQD_EVAL(z);
        }

        factor = tmpl_Float_Exp(-abs_x*abs_x - 0.5625F + P / Q);
        erf_abs_x = 1.0F - factor * rcpr_abs_x;
    }

    /*  The error function is odd. If x is negative, return minus the result. */
    if (x < 0.0F)
        return -erf_abs_x;

    /*  Otherwise we have computed the error function for x and can return.   */
    return erf_abs_x;
}
/*  End of Erf_Float_Sun.                                                     */

/*  Undefine everything in case some wants to #include this file.             */
#undef PA0
#undef PA1
#undef PA2

#undef QA0
#undef QA1
#undef QA2
#undef QA3

#undef PB0
#undef PB1
#undef PB2
#undef PB3

#undef QB0
#undef QB1
#undef QB2
#undef QB3
#undef QB4

#undef PC0
#undef PC1
#undef PC2
#undef PC3

#undef QC0
#undef QC1
#undef QC2
#undef QC3
#undef QC4

#undef PD0
#undef PD1
#undef PD2
#undef PD3
#undef PD4

#undef QD0
#undef QD1
#undef QD2
#undef QD3
#undef QD4

#undef POLYPA_EVAL
#undef POLYPB_EVAL
#undef POLYPC_EVAL
#undef POLYPD_EVAL

#undef POLYQA_EVAL
#undef POLYQB_EVAL
#undef POLYQC_EVAL
#undef POLYQD_EVAL
