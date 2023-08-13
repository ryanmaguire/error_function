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
#define C0 (+8.45062911510467529297E-01)

/*  Scale factor for the leading term of the Taylor series.                   */
#define TWO_BY_ROOT_PI (+1.1283791670955125738961589031215451716881012586E+00)

/*  Coefficients for the rational approximation to erf on [0, 0.84375].       */
#define PA0 (+1.28379167095512558561E-01)
#define PA1 (-3.25042107247001499370E-01)
#define PA2 (-2.84817495755985104766E-02)
#define PA3 (-5.77027029648944159157E-03)
#define PA4 (-2.37630166566501626084E-05)

#define QA0 (+1.00000000000000000000E+00)
#define QA1 (+3.97917223959155352819E-01)
#define QA2 (+6.50222499887672944485E-02)
#define QA3 (+5.08130628187576562776E-03)
#define QA4 (+1.32494738004321644526E-04)
#define QA5 (-3.96022827877536812320E-06)

/*  Helper macros for evaluating polynomials using Horner's method.           */
#define POLYPA_EVAL(t) \
PA0 + t*(PA1 + t*(PA2 + t*(PA3 + t*PA4)))

#define POLYQA_EVAL(t) \
QA0 + t*(QA1 + t*(QA2 + t*(QA3 + t*(QA4 + t*QA5))))

/*  Coefficients for the rational approximation to erf on [0.84375, 1.25].    */
#define PB0 (-2.36211856075265944077E-03)
#define PB1 (+4.14856118683748331666E-01)
#define PB2 (-3.72207876035701323847E-01)
#define PB3 (+3.18346619901161753674E-01)
#define PB4 (-1.10894694282396677476E-01)
#define PB5 (+3.54783043256182359371E-02)
#define PB6 (-2.16637559486879084300E-03)

#define QB0 (+1.00000000000000000000E+00)
#define QB1 (+1.06420880400844228286E-01)
#define QB2 (+5.40397917702171048937E-01)
#define QB3 (+7.18286544141962662868E-02)
#define QB4 (+1.26171219808761642112E-01)
#define QB5 (+1.36370839120290507362E-02)
#define QB6 (+1.19844998467991074170E-02)

/*  Helper macros for evaluating polynomials using Horner's method.           */
#define POLYPB_EVAL(t) \
PB0 + t*(PB1 + t*(PB2 + t*(PB3 + t*(PB4 + t*(PB5 + t*PB6)))))

#define POLYQB_EVAL(t) \
QB0 + t*(QB1 + t*(QB2 + t*(QB3 + t*(QB4 + t*(QB5 + t*QB6)))))

/*  Coefficients for the exponential approximation to erf on [1.25, 2.8571].  */
#define PC0 (-9.86494403484714822705E-03)
#define PC1 (-6.93858572707181764372E-01)
#define PC2 (-1.05586262253232909814E+01)
#define PC3 (-6.23753324503260060396E+01)
#define PC4 (-1.62396669462573470355E+02)
#define PC5 (-1.84605092906711035994E+02)
#define PC6 (-8.12874355063065934246E+01)
#define PC7 (-9.81432934416914548592E+00)

#define QC0 (+1.00000000000000000000E+00)
#define QC1 (+1.96512716674392571292E+01)
#define QC2 (+1.37657754143519042600E+02)
#define QC3 (+4.34565877475229228821E+02)
#define QC4 (+6.45387271733267880336E+02)
#define QC5 (+4.29008140027567833386E+02)
#define QC6 (+1.08635005541779435134E+02)
#define QC7 (+6.57024977031928170135E+00)
#define QC8 (-6.04244152148580987438E-02)

/*  Helper macros for evaluating polynomials using Horner's method.           */
#define POLYPC_EVAL(t) \
PC0 + t*(PC1 + t*(PC2 + t*(PC3 + t*(PC4 + t*(PC5 + t*(PC6 + t*PC7))))))

#define POLYQC_EVAL(t) \
QC0+t*(QC1+t*(QC2+t*(QC3+t*(QC4+t*(QC5+t*(QC6+t*(QC7+t*QC8)))))))

/*  Coefficients for the exponential approximation to erf on [2.8571, 6.0].   */
#define PD0 (-9.86494292470009928597E-03)
#define PD1 (-7.99283237680523006574E-01)
#define PD2 (-1.77579549177547519889E+01)
#define PD3 (-1.60636384855821916062E+02)
#define PD4 (-6.37566443368389627722E+02)
#define PD5 (-1.02509513161107724954e+03)
#define PD6 (-4.83519191608651397019e+02)

#define QD0 (+1.00000000000000000000E+00)
#define QD1 (+3.03380607434824582924E+01)
#define QD2 (+3.25792512996573918826E+02)
#define QD3 (+1.53672958608443695994E+03)
#define QD4 (+3.19985821950859553908e+03)
#define QD5 (+2.55305040643316442583E+03)
#define QD6 (+4.74528541206955367215E+02)
#define QD7 (-2.24409524465858183362e+01)

/*  Helper macros for evaluating polynomials using Horner's method.           */
#define POLYPD_EVAL(t) \
PD0 + t*(PD1 + t*(PD2 + t*(PD3 + t*(PD4 + t*(PD5 + t*PD6)))))

#define POLYQD_EVAL(t) \
QD0 + t*(QD1 + t*(QD2 + t*(QD3 + t*(QD4 + t*(QD5 + t*(QD6 + t*QD7))))))

/*  Computes erf(x) using the algorithm outlined in the comments of SunOS's   *
 *  implementation of libm. This does not use their actual code, which makes  *
 *  use of macros assuming IEEE-754 support. Instead this uses the same idea  *
 *  but with code that is compliant with the C89 standard and quite portable. */
double Erf_Double_Sun(double x)
{
    /*  erf(x) is odd. We'll compute erf(|x|) and check the sign of x later.  */
    double erf_abs_x;
    const double abs_x = tmpl_Double_Abs(x);

    /*  For |x| >= 6 we have erf(|x|) = 1 to double precision. Return this.   */
    if (abs_x >= 6.0)
        erf_abs_x = 1.0;

    /*  For very small values, use the leading term of the Taylor series.     */
    else if (abs_x < 3.725290298461914e-09)
        return TWO_BY_ROOT_PI*x;

    /*  Slightly larger, use a Remez rational minimax approximation.          */
    else if (abs_x < 0.84375)
    {
        const double z = abs_x * abs_x;
        const double P = POLYPA_EVAL(z);
        const double Q = POLYQA_EVAL(z);
        return x + x * P/Q;
    }

    /*  For larger values, use a rational approximation centered at x = 1.    */
    else if (abs_x < 1.25)
    {
        const double z = abs_x - 1.0;
        const double P = POLYPB_EVAL(z);
        const double Q = POLYQB_EVAL(z);
        erf_abs_x = C0 + P/Q;
    }

    /*  And for all other values, use an exponentiated rational approximation.*/
    else
    {
        const double rcpr_abs_x = 1.0 / abs_x;
        const double z = rcpr_abs_x * rcpr_abs_x;
        double P, Q, factor;

        if (abs_x < 2.857142857142857)
        {
            P = POLYPC_EVAL(z);
            Q = POLYQC_EVAL(z);
        }

        else
        {
            P = POLYPD_EVAL(z);
            Q = POLYQD_EVAL(z);
        }

        factor = tmpl_Double_Exp(-abs_x*abs_x - 0.5625 + P / Q);
        erf_abs_x = 1.0 - factor * rcpr_abs_x;
    }

    /*  The error function is odd. If x is negative, return minus the result. */
    if (x < 0.0)
        return -erf_abs_x;

    /*  Otherwise we have computed the error function for x and can return.   */
    return erf_abs_x;
}
/*  End of Erf_Double_Sun.                                                    */

/*  Undefine everything in case some wants to #include this file.             */
#undef PA0
#undef PA1
#undef PA2
#undef PA3
#undef PA4

#undef QA0
#undef QA1
#undef QA2
#undef QA3
#undef QA4
#undef QA5

#undef PB0
#undef PB1
#undef PB2
#undef PB3
#undef PB4
#undef PB5
#undef PB6

#undef QB0
#undef QB1
#undef QB2
#undef QB3
#undef QB4
#undef QB5
#undef QB6

#undef PC0
#undef PC1
#undef PC2
#undef PC3
#undef PC4
#undef PC5
#undef PC6
#undef PC7

#undef QC0
#undef QC1
#undef QC2
#undef QC3
#undef QC4
#undef QC5
#undef QC6
#undef QC7
#undef QC8

#undef PD0
#undef PD1
#undef PD2
#undef PD3
#undef PD4
#undef PD5
#undef PD6

#undef QD0
#undef QD1
#undef QD2
#undef QD3
#undef QD4
#undef QD5
#undef QD6
#undef QD7

#undef POLYPA_EVAL
#undef POLYPB_EVAL
#undef POLYPC_EVAL
#undef POLYPD_EVAL

#undef POLYQA_EVAL
#undef POLYQB_EVAL
#undef POLYQC_EVAL
#undef POLYQD_EVAL
