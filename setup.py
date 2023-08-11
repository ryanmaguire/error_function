"""
################################################################################
#                                   LICENSE                                    #
################################################################################
#   This file is part of error_function.                                       #
#                                                                              #
#   error_function is free software: you can redistribute it and/or modify it  #
#   under the terms of the GNU General Public License as published by          #
#   the Free Software Foundation, either version 3 of the License, or          #
#   (at your option) any later version.                                        #
#                                                                              #
#   error_function is distributed in the hope that it will be useful,          #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with error_function.  If not, see <https://www.gnu.org/licenses/>.   #
################################################################################
"""

import distutils.core
import numpy

srclist = [
    "src/setup.c",
    "src/erf_abramowitz_and_stegun_double.c",
    "src/erf_abramowitz_and_stegun_float.c",
    "src/erf_abramowitz_and_stegun_ldouble.c",
    "src/erf_abramowitz_and_stegun_rational_double.c",
    "src/erf_abramowitz_and_stegun_rational_float.c",
    "src/erf_abramowitz_and_stegun_rational_ldouble.c",
    "src/erf_winitzki_double.c",
    "src/erf_winitzki_float.c",
    "src/erf_winitzki_ldouble.c"
]

distutils.core.setup(
    name = 'erfpy',
    version = '0.1',
    description = 'Various algorithms for the real-valued error function',
    author='Ryan Maguire',
    ext_modules = [
        distutils.core.Extension(
            'erfpy',
            srclist,
            include_dirs = [numpy.get_include()],
            libraries = ['tmpl']
        )
    ]
)
