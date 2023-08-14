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
#define C0 (+8.45062911510467529296875E-01L)

/*  Scale factor for the leading term of the Taylor series.                   */
#define TWO_BY_ROOT_PI (+1.1283791670955125738961589031215451716881012586E+00L)

/*  Coefficients for the rational approximation to erf on [0, 0.84375].       */
#define PA0 (+1.283791670955125738886075e-01L)
#define PA1 (-3.211371473841651581477945e-01L)
#define PA2 (-3.789711267118296883266073e-02L)
#define PA3 (-7.831016249530390647944696e-03L)
#define PA4 (-3.039038177926934657287030e-04L)
#define PA5 (-1.892737494971971821647348e-05L)

#define QA0 (+1.000000000000000000000000e+00L)
#define QA1 (+4.283346191735379040121338e-01L)
#define QA2 (+8.079912388447361625226726e-02L)
#define QA3 (+8.517355594538185460097254e-03L)
#define QA4 (+5.158696927867366838085203e-04L)
#define QA5 (+1.565249677548820031842715e-05L)
#define QA6 (+1.143433646151683005684191e-07L)

/*  Helper macros for evaluating polynomials using Horner's method.           */
#define POLYPA_EVAL(t) \
PA0+t*(PA1+t*(PA2+t*(PA3+t*(PA4+t*PA5))))

#define POLYQA_EVAL(t) \
QA0+t*(QA1+t*(QA2+t*(QA3+t*(QA4+t*(QA5+t*QA6)))))

/*  Coefficients for the rational approximation to erf on [0.84375, 1.25].    */
#define PB0 (-2.362118560752659955555513e-03L)
#define PB1 (+4.134033634759403336243061e-01L)
#define PB2 (-1.171056183678228933866790e-01L)
#define PB3 (+9.729445257088016520835572e-02L)
#define PB4 (+3.691866973611315094370641e-02L)
#define PB5 (-5.176793364366149594906803e-03L)
#define PB6 (+4.062563959186227983980474e-03L)
#define PB7 (+2.060638519506302803470926e-04L)

#define QB0 (+9.999999999999999999457899e-01L)
#define QB1 (+7.214430185550756159901069e-01L)
#define QB2 (+6.242391621170184547043661e-01L)
#define QB3 (+3.067854689398429930040704e-01L)
#define QB4 (+1.329203815108533268097256e-01L)
#define QB5 (+4.559279314400553275807445e-02L)
#define QB6 (+1.017986985805487830452850e-02L)
#define QB7 (+2.193336601938737913289240e-03L)

/*  Helper macros for evaluating polynomials using Horner's method.           */
#define POLYPB_EVAL(t) \
PB0+t*(PB1+t*(PB2+t*(PB3+t*(PB4+t*(PB5+t*(PB6+t*PB7))))))

#define POLYQB_EVAL(t) \
QB0+t*(QB1+t*(QB2+t*(QB3+t*(QB4+t*(QB5+t*(QB6+t*QB7))))))

/*  Coefficients for the exponential approximation to erf on [1.25, 2.8571].  */
#define PC0 (-9.864943092668513353233099e-03L)
#define PC1 (-7.366355527888176743566209e-01L)
#define PC2 (-1.347354008051589297108308e+01L)
#define PC3 (-1.021261197225424459528176e+02L)
#define PC4 (-3.681385465941502179443301e+02L)
#define PC5 (-6.459288018931153199098105e+02L)
#define PC6 (-5.255573741023815239059580e+02L)
#define PC7 (-1.727270042388650067161082e+02L)
#define PC8 (-1.606758080227204993004064e+01L)

#define QC0 (+9.999999999999999999457899e-01L)
#define QC1 (+2.398752552695060142294581e+01L)
#define QC2 (+2.133592090764565437771560e+02L)
#define QC3 (+9.018889241674588011510316e+02L)
#define QC4 (+1.933882820697022772105100e+03L)
#define QC5 (+2.083777771258359691630702e+03L)
#define QC6 (+1.049460201650984270460576e+03L)
#define QC7 (+2.079632269886748182907699e+02L)
#define QC8 (+1.014474124618066440232933e+01L)
#define QC9 (-7.234661769911255643349546e-02L)

/*  Helper macros for evaluating polynomials using Horner's method.           */
#define POLYPC_EVAL(t) \
PC0+t*(PC1+t*(PC2+t*(PC3+t*(PC4+t*(PC5+t*(PC6+t*(PC7+t*PC8)))))))

#define POLYQC_EVAL(t) \
QC0+t*(QC1+t*(QC2+t*(QC3+t*(QC4+t*(QC5+t*(QC6+t*(QC7+t*(QC8+t*QC9))))))))

/*  Coefficients for the exponential approximation to erf on [2.8571, 6.6].   */
#define PD0 (-9.864942924704175847177884e-03L)
#define PD1 (-8.164487895849486653182536e-01L)
#define PD2 (-1.911251735369680328918696e+01L)
#define PD3 (-1.887875083867998777520381e+02L)
#define PD4 (-8.657956700744415108883523e+02L)
#define PD5 (-1.791295223806467901916939e+03L)
#define PD6 (-1.432100914726251923436706e+03L)
#define PD7 (-2.838647383144148992339950e+02L)

#define QD0 (+1.000000000000000000000000e+00L)
#define QD1 (+3.207811661981680352237478e+01L)
#define QD2 (+3.749093199818248671184584e+02L)
#define QD3 (+2.011162638369640094926361e+03L)
#define QD4 (+5.128720611187760284988002e+03L)
#define QD5 (+5.813623703138272055479518e+03L)
#define QD6 (+2.394648385194624399918695e+03L)
#define QD7 (+2.025827286619646666193661e+02L)

/*  Helper macros for evaluating polynomials using Horner's method.           */
#define POLYPD_EVAL(t) \
PD0+t*(PD1+t*(PD2+t*(PD3+t*(PD4+t*(PD5+t*(PD6+t*PD7))))))

#define POLYQD_EVAL(t) \
QD0+t*(QD1+t*(QD2+t*(QD3+t*(QD4+t*(QD5+t*(QD6+t*QD7))))))

/*  Computes erf(x) using the algorithm outlined in the comments of SunOS's   *
 *  implementation of libm. This does not use their actual code, which makes  *
 *  use of macros assuming IEEE-754 support. Instead this uses the same idea  *
 *  but with code that is compliant with the C89 standard and quite portable. */
long double Erf_LDouble_Sun(long double x)
{
    /*  erf(x) is odd. We'll compute erf(|x|) and check the sign of x later.  */
    long double erf_abs_x;
    const long double abs_x = tmpl_LDouble_Abs(x);

    /*  For |x| >= 6.6 we have erf(|x|) = 1 to extended precision.            */
    if (abs_x >= 6.6666259765625L)
        erf_abs_x = 1.0L;

    /*  For very small values, use the leading term of the Taylor series.     */
    else if (abs_x < 1.1641532182693481E-10L)
        return TWO_BY_ROOT_PI*x;

    /*  Slightly larger, use a Remez rational minimax approximation.          */
    else if (abs_x < 0.84375L)
    {
        const long double z = abs_x * abs_x;
        const long double P = POLYPA_EVAL(z);
        const long double Q = POLYQA_EVAL(z);
        return x + x * P/Q;
    }

    /*  For larger values, use a rational approximation centered at x = 1.    */
    else if (abs_x < 1.25L)
    {
        const long double z = abs_x - 1.0L;
        const long double P = POLYPB_EVAL(z);
        const long double Q = POLYQB_EVAL(z);
        erf_abs_x = C0 + P/Q;
    }

    /*  And for all other values, use an exponentiated rational approximation.*/
    else
    {
        const long double rcpr_abs_x = 1.0L / abs_x;
        const long double z = rcpr_abs_x * rcpr_abs_x;
        long double P, Q, factor;

        if (abs_x < 2.857142857142857L)
        {
            P = POLYPC_EVAL(z);
            Q = POLYQC_EVAL(z);
        }

        else
        {
            P = POLYPD_EVAL(z);
            Q = POLYQD_EVAL(z);
        }

        factor = tmpl_LDouble_Exp(-abs_x*abs_x - 0.5625L + P / Q);
        erf_abs_x = 1.0L - factor * rcpr_abs_x;
    }

    /*  The error function is odd. If x is negative, return minus the result. */
    if (x < 0.0L)
        return -erf_abs_x;

    /*  Otherwise we have computed the error function for x and can return.   */
    return erf_abs_x;
}
/*  End of Erf_LDouble_Sun.                                                   */

/*  Undefine everything in case some wants to #include this file.             */
#undef PA0
#undef PA1
#undef PA2
#undef PA3
#undef PA4
#undef PA5

#undef QA0
#undef QA1
#undef QA2
#undef QA3
#undef QA4
#undef QA5
#undef QA6

#undef PB0
#undef PB1
#undef PB2
#undef PB3
#undef PB4
#undef PB5
#undef PB6
#undef PB7

#undef QB0
#undef QB1
#undef QB2
#undef QB3
#undef QB4
#undef QB5
#undef QB6
#undef QB7

#undef PC0
#undef PC1
#undef PC2
#undef PC3
#undef PC4
#undef PC5
#undef PC6
#undef PC7
#undef PC8

#undef QC0
#undef QC1
#undef QC2
#undef QC3
#undef QC4
#undef QC5
#undef QC6
#undef QC7
#undef QC8
#undef QC8
#undef QC9

#undef PD0
#undef PD1
#undef PD2
#undef PD3
#undef PD4
#undef PD5
#undef PD6
#undef PD7

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
