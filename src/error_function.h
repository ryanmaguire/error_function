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

/*  Include guard to prevent including this file twice.                       */
#ifndef ERROR_FUNCTION_H
#define ERROR_FUNCTION_H

extern double Erf_Double_Abramowitz_and_Stegun(double x);
extern float Erf_Float_Abramowitz_and_Stegun(float x);
extern long double Erf_LDouble_Abramowitz_and_Stegun(long double x);

extern double Erf_Double_Abramowitz_and_Stegun_Rational(double x);
extern float Erf_Float_Abramowitz_and_Stegun_Rational(float x);
extern long double Erf_LDouble_Abramowitz_and_Stegun_Rational(long double x);

extern double Erf_Double_Karagiannidis_and_Lioumpas(double x);
extern float Erf_Float_Karagiannidis_and_Lioumpas(float x);
extern long double Erf_LDouble_Karagiannidis_and_Lioumpas(long double x);

extern double Erf_Double_Numerical_Recipes(double x);
extern float Erf_Float_Numerical_Recipes(float x);
extern long double Erf_LDouble_Numerical_Recipes(long double x);

extern double Erf_Double_Sun(double x);

extern double Erf_Double_Winitzki(double x);
extern float Erf_Float_Winitzki(float x);
extern long double Erf_LDouble_Winitzki(long double x);

#endif
/*  End of include guard.                                                     */
