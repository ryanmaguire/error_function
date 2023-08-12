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

import os
import setuptools
import numpy

srclist = []

for file in os.listdir("src/"):
    if file[-1] == "c":
        srclist.append("src/" + file)

setuptools.setup(
    name = 'erfpy',
    version = '0.1',
    description = 'Various algorithms for the real-valued error function',
    author='Ryan Maguire',
    ext_modules = [
        setuptools.Extension(
            'erfpy',
            srclist,
            include_dirs = [numpy.get_include()],
            libraries = ['tmpl']
        )
    ]
)
