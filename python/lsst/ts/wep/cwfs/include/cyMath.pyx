# This file is part of ts_wep.
#
# Developed for the LSST Telescope and Site Systems.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# -*- coding: utf-8 -*-

# This is to optimize the speed of calculation related to annular Zernike polynomials

cdef extern from "math.h":
    double sqrt(double theta)
    double sin(double theta)
    double cos(double theta)
    double atan2(double theta, double theta)

from numpy cimport ndarray
import numpy as np
cimport numpy as np
cimport cython


@cython.boundscheck(False)
def ZernikeAnnularJacobian(ndarray[np.float64_t, ndim=1] Z not None,
                           ndarray[np.float64_t, ndim=1] x not None,
                           ndarray[np.float64_t, ndim=1] y not None, e, atype):

    cdef Py_ssize_t n = x.shape[0]
    cdef ndarray[np.float64_t, ndim=1] out = np.zeros_like(x)

    _ZernikeAnnularJacobian(&out[0], &Z[0], &x[0], &y[0], e, atype, n)

    return out


cdef _ZernikeAnnularJacobian(double *out, double *Z, double *x, double *y, double e,
                             str atype, int n):

    # Parameters of constant
    cdef double e2 = e**2
    cdef double e4 = e2*e2
    cdef double e6 = e4*e2
    cdef double e8 = e6*e2
    cdef double e10 = e8*e2
    cdef double e12 = e10*e2
    cdef double e14 = e12 * e2
    cdef double e16 = e14 * e2

    cdef double sqrt_3 = sqrt(3)
    cdef double sqrt_5 = sqrt(5)
    cdef double sqrt_6 = sqrt(6)
    cdef double sqrt_7 = sqrt(7)
    cdef double sqrt_8 = sqrt(8)
    cdef double sqrt_10 = sqrt(10)
    cdef double sqrt_12 = sqrt(12)

    # 1st order
    cdef double den1 = 1 - e2
    cdef double den2 = sqrt((1 - e2)**2 * (1 + e2) * (1 + 4 * e2 + e4))
    cdef double den3 = (1 - e2)**2

    cdef double den4 = (1 - e2)**3 * (1 + e2 + e4)
    cdef double num4 = sqrt((1 - e2)**4 * (1 + e2 + e4) /
                            (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8))

    cdef double den5 = (1 - e2)**3 * (1 + 4 * e2 + e4)
    cdef double num5 = sqrt((1 - e2)**2 * (1 + 4 * e2 + e4) /
                            (1 + 9 * e2 + 9 * e4 + e6))

    cdef double den6 = (1 - e2)**4 * (1 + e2) * (1 + e4)
    cdef double num6 = sqrt((1 - e2)**6 * (1 + e2) * (1 + e4) /
                            (1 + 4 * e2 + 10 * e4 + 20 * e6 +
                            10 * e8 + 4 * e10 + e12))

    cdef double den7 = (1 - e2)**3

    # 2nd order
    cdef double den2_2 = (1 + e2 + e4)
    cdef double den2_3 = (1 - e2)**2 * (1 + e2) * (1 + 4 * e2 + e4)
    cdef double den2_4 = (1 + e2 + e4 + e6)
    cdef double den2_5 = (1 - e2)**4

    cdef double den2_6 = (1 - e2)**6 * (1 + e2 + e4)**2
    cdef double num2_6 = ((1 - e2)**4 * (1 + e2 + e4) /
                          (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8))

    cdef double den2_7 = (1 + e2 + e4 + e6 + e8)

    cdef double den2_8 = (1 - e2)**6 * (1 + 4 * e2 + e4)**2
    cdef double num2_8 = (1 - e2)**2 * (1 + 4 * e2 + e4) / (1 + 9 * e2 + 9 * e4 + e6)

    cdef double den2_9 = (1 - e2)**8 * (1 + e2)**2 * (1 + e4)**2
    cdef double num2_9 = (1 - e2)**6 * (1 + e2) * (1 + e4) / \
                         (1 + 4 * e2 + 10 * e4 + 20 * e6 + 10 * e8 + 4 * e10 + e12)

    cdef double den2_10 = (1 + e2 + e4 + e6 + e8 + e10)
    cdef double den2_11 = (1 - e2)**6

    # Parameters in loop
    cdef double x2, y2, x4, y4, xy, r2, r4, x6, y6, temp, x_c, y_c
    cdef int ii

    if (atype == "1st"):

        for ii in range(n):

            x_c = x[ii]
            y_c = y[ii]

            x2 = x_c * x_c
            y2 = y_c * y_c
            xy = x_c * y_c
            r2 = x2 + y2
            x4 = x2 * x2
            y4 = y2 * y2

            temp = Z[0] * 0 * x_c  # to make d an array with the same size as x
            temp += Z[1] * 0
            temp += Z[2] * 0

            temp += Z[3] * sqrt_3 * 8 / den1
            temp += Z[4] * sqrt_6 * 0
            temp += Z[5] * sqrt_6 * 0

            temp += Z[6] * sqrt_8 * 24 * y_c * (1 + e2) / den2
            temp += Z[7] * sqrt_8 * 24 * x_c * (1 + e2) / den2
            temp += Z[8] * sqrt_8 * 0
            temp += Z[9] * sqrt_8 * 0

            temp += Z[10] * sqrt_5 * (96 * r2 - 24 * (1 + e2)) / den3

            temp += Z[11] * sqrt_10 * 48 * (x2 - y2) * \
                (1 + e2 + e4) * num4 / den4
            temp += Z[12] * sqrt_10 * 96 * xy * (1 + e2 + e4) * num4 / den4
            temp += Z[13] * sqrt_10 * 0
            temp += Z[14] * sqrt_10 * 0

            temp += Z[15] * sqrt_12 * 48 * x_c * (
                5 * r2 * (1 + 4 * e2 + e4) - 2 *
                (1 + 4 * e2 + 4 * e4 + e6)) * num5 / den5
            temp += Z[16] * sqrt_12 * 48 * y_c * (
                5 * r2 * (1 + 4 * e2 + e4) - 2 *
                (1 + 4 * e2 + 4 * e4 + e6)) * num5 / den5

            temp += Z[17] * sqrt_12 * 80.0 * x_c * \
                (x2 - 3.0 * y2) * (1 + e2) * (1 + e4) * num6 / den6
            temp += Z[18] * sqrt_12 * 80.0 * y_c * \
                (3 * x2 - y2) * (1 + e2) * (1 + e4) * num6 / den6
            temp += Z[19] * sqrt_12 * 0
            temp += Z[20] * sqrt_12 * 0

            temp += Z[21] * sqrt_7 * 48 * (
                e4 - 10 * e2 * x2 - 10 * e2 * y2 +
                3 * e2 + 15 * x4 + 30 * x2 * y2 - 10 * x2 +
                15 * y4 - 10 * y2 + 1) / den7

            out[ii] = temp

    elif (atype == "2nd"):

        for ii in range(n):

            x_c = x[ii]
            y_c = y[ii]

            x2 = x_c * x_c
            y2 = y_c * y_c
            xy = x_c * y_c
            r2 = x2 + y2
            x4 = x2 * x2
            x6 = x4 * x2
            y4 = y2 * y2
            y6 = y4 * y2

            temp = Z[0]**2 * 0 * x_c  # to make d an array with the same size as x
            temp += Z[1]**2 * 0
            temp += Z[2]**2 * 0

            temp += Z[3]**2 * (3) * 16 / den1 / den1

            temp += Z[4]**2 * (6) * (-4) / den2_2
            temp += Z[5]**2 * (6) * (-4) / den2_2

            temp += Z[6]**2 * (8) * (108 * y2 - 36 * x2) * (1 + e2) / den2_3
            temp += Z[7]**2 * (8) * (108 * x2 - 36 * y2) * (1 + e2) / den2_3

            temp += Z[8]**2 * (8) * (-36 * r2) / den2_4
            temp += Z[9]**2 * (8) * (-36 * r2) / den2_4

            temp += Z[10]**2 * (5) * 144 * (1 + e2 - 2 * r2) * \
                (1 + e2 - 6 * r2) / den2_5

            temp += Z[11]**2 * (10) * 36 * (
                8 * (1 + e2 + e4) * x2 - 1 - e2 - e4 - e6) * \
                (1 + e2 + e4 + e6 - 8 * (1 + e2 + e4) * y2) * num2_6 / den2_6
            temp += Z[12]**2 * (10) * 36 * (
                -4 * (x_c - y_c)**2 * (e4 + e2 + 1) + 1 + e2 + e4 + e6) * \
                (4 * (x_c + y_c)**2 * (e4 + e2 + 1) - 1 - e2 - e4 - e6) * num2_6 / den2_6

            temp += Z[13]**2 * (10) * (-144) * r2**2 / den2_7
            temp += Z[14]**2 * (10) * (-144) * r2**2 / den2_7

            temp += Z[15]**2 * (12) * 64 * (
                (3 * e6 - 5 * e4 * r2 + 12 * e4 - 20 * e2 * r2 +
                 12 * e2 - 5 * r2 + 3) *
                (9 * e6 * x2 - 3 * e6 * y2 - 25 * e4 * x4 - 20 * e4 * x2 * y2 +
                 36 * e4 * x2 + 5 * e4 * y4 - 12 * e4 * y2 - 100 * e2 * x4 -
                 80 * e2 * x2 * y2 + 36 * e2 * x2 + 20 * e2 * y4 -
                 12 * e2 * y2 - 25 * x4 - 20 * x2 * y2 +
                 9 * x2 + 5 * y4 - 3 * y2)) * num2_8 / den2_8
            temp += Z[16]**2 * (12) * 64 * (
                -(3 * e6 - 5 * e4 * r2 + 12 * e4 - 20 * e2 * r2 + 12 * e2 -
                  5 * r2 + 3) * (3 * e6 * x2 - 9 * e6 * y2 - 5 * e4 * x4 +
                                 20 * e4 * x2 * y2 + 12 * e4 * x2 + 25 * e4 * y4 -
                                 36 * e4 * y2 - 20 * e2 * x4 + 80 * e2 * x2 * y2 +
                                 12 * e2 * x2 + 100 * e2 * y4 -
                                 36 * e2 * y2 - 5 * x4 +
                                 20 * x2 * y2 + 3 * x2 + 25 * y4 -
                                 9 * y2)) * num2_8 / den2_8

            temp += Z[17]**2 * (12) * 16.0 * (
                - 36 * e16 * x2 - 36 * e16 * y2 + 180 * e14 * x4 +
                360 * e14 * x2 * y2 - 72 * e14 * x2 +
                180 * e14 * y4 - 72 * e14 * y2 - 125 * e12 * x6 -
                1275 * e12 * x4 * y2 + 360 * e12 * x4 + 225 * e12 * x2 * y4 +
                720 * e12 * x2 * y2 - 108 * e12 * x2 - 225 * e12 *
                y6 + 360 * e12 * y4 - 108 * e12 * y2 - 250 * e10 * x6 -
                2550 * e10 * x4 * y2 + 540 * e10 * x4 + 450 * e10 * x2 *
                y4 + 1080 * e10 * x2 * y2 - 144 * e10 * x2 - 450 * e10 * y6 +
                540 * e10 * y4 - 144 * e10 * y2 - 375 * e8 * x6 - 3825 *
                e8 * x4 * y2 + 720 * e8 * x4 + 675 * e8 * x2 * y4 +
                1440 * e8 * x2 * y2 - 180 * e8 * x2 - 675 * e8 * y6 + 720 *
                e8 * y4 - 180 * e8 * y2 - 500 * e6 * x6 - 5100 * e6 * x4 * y2 +
                720 * e6 * x4 + 900 * e6 * x2 * y4 + 1440 * e6 * x2 * y2 -
                144 * e6 * x2 - 900 * e6 * y6 + 720 * e6 * y4 - 144 * e6 * y2 -
                375 * e4 * x6 - 3825 * e4 * x4 * y2 + 540 * e4 * x4 + 675 * e4 *
                x2 * y4 + 1080 * e4 * x2 * y2 - 108 * e4 * x2 - 675 * e4 * y6 +
                540 * e4 * y4 - 108 * e4 * y2 - 250 * e2 * x6 - 2550 * e2 * x4 *
                y2 + 360 * e2 * x4 + 450 * e2 * x2 * y4 + 720 * e2 * x2 * y2 -
                72 * e2 * x2 - 450 * e2 * y6 + 360 * e2 * y4 - 72 * e2 *
                y2 - 125 * x6 - 1275 * x4 * y2 + 180 * x4 + 225 * x2 * y4 +
                360 * x2 * y2 - 36 * x2 - 225 * y6 + 180 * y4 -
                36 * y2) * num2_9 / den2_9
            temp += Z[18]**2 * (12) * 16.0 * ((
                - 225 * e12 - 450 * e10 - 675 * e8 - 900 * e6 - 675 * e4 -
                450 * e2 - 225) * x6 +
                (180 * e14 + 225 * e12 * y2 + 360 * e12 + 450 * e10 * y2 +
                 540 * e10 + 675 * e8 * y2 + 720 * e8 + 900 * e6 * y2 +
                 720 * e6 + 675 * e4 * y2 + 540 * e4 + 450 * e2 * y2 +
                 360 * e2 + 225 * y2 + 180) * x4 +
                (- 36 * e16 + 360 * e14 * y2 - 72 * e14 - 1275 * e12 * y4 +
                 720 * e12 * y2 - 108 * e12 - 2550 * e10 * y4 +
                 1080 * e10 * y2 - 144 * e10 - 3825 * e8 * y4 + 1440 *
                 e8 * y2 - 180 * e8 - 5100 * e6 * y4 + 1440 * e6 * y2 -
                 144 * e6 - 3825 * e4 * y4 + 1080 * e4 * y2 - 108 * e4 -
                 2550 * e2 * y4 + 720 * e2 * y2 - 72 * e2 - 1275 * y4 +
                 360 * y2 - 36) * x2 - 36 * e16 * y2 + 180 * e14 * y4 -
                72 * e14 * y2 - 125 * e12 * y6 + 360 * e12 * y4 - 108 * e12 * y2 -
                250 * e10 * y6 + 540 * e10 * y4 - 144 * e10 * y2 - 375 *
                e8 * y6 + 720 * e8 * y4 - 180 * e8 * y2 - 500 * e6 * y6 +
                720 * e6 * y4 - 144 * e6 * y2 - 375 * e4 * y6 + 540 * e4 * y4 -
                108 * e4 * y2 - 250 * e2 * y6 + 360 * e2 * y4 - 72 * e2 * y2 -
                125 * y6 + 180 * y4 - 36 * y2) * num2_9 / den2_9

            temp += Z[19]**2 * (12) * (-400) * r2**3 / den2_10
            temp += Z[20]**2 * (12) * (-400) * r2**3 / den2_10

            temp += Z[21]**2 * (7) * 576 * ((
                e4 - 5 * e2 * x2 - 5 * e2 * y2 + 3 * e2 + 5 * x4 +
                10 * x2 * y2 - 5 * x2 + 5 * y4 - 5 * y2 + 1) *
                (e4 - 15 * e2 * x2 - 15 * e2 * y2 + 3 * e2 +
                 25 * x4 + 50 * x2 * y2 - 15 * x2 + 25 * y4 - 15 * y2 + 1)) / den2_11

            out[ii] = temp


def ZernikeAnnularGrad(ndarray[np.float64_t, ndim=1] Z not None, ndarray[np.float64_t, ndim=1] x not None, 
                       ndarray[np.float64_t, ndim=1] y not None, e, axis):

    cdef Py_ssize_t n = x.shape[0]
    cdef ndarray[np.float64_t, ndim=1] d = np.zeros_like(x)

    _ZernikeAnnularGrad(&d[0], &Z[0], &x[0], &y[0], e, axis, n)

    return d


cdef _ZernikeAnnularGrad(double *d, double *Z, double *x, double *y, double e,
                         str axis, int n):

    # Parameters of constant
    cdef double e2 = e**2
    cdef double e4 = e2*e2
    cdef double e6 = e4*e2
    cdef double e8 = e6*e2
    cdef double e10 = e8*e2
    cdef double e12 = e10*e2

    cdef double sqrt_3 = sqrt(3)
    cdef double sqrt_5 = sqrt(5)
    cdef double sqrt_6 = sqrt(6)
    cdef double sqrt_7 = sqrt(7)
    cdef double sqrt_8 = sqrt(8)
    cdef double sqrt_10 = sqrt(10)
    cdef double sqrt_12 = sqrt(12)

    cdef double den1 = sqrt(1 + e2)
    cdef double den2 = 1 - e2
    cdef double den3 = sqrt(1 + e2 + e4)
    cdef double den4 = sqrt((1 - e2)**2 * (1 + e2) * (1 + 4 * e2 + e4))
    cdef double den5 = sqrt(1 + e2 + e4 + e6)
    cdef double den6 = (1 - e2)**2

    cdef double den7 = (1 - e2)**3 * (1 + e2 + e4)
    cdef double num7 = sqrt((1 - e2)**4 * (1 + e2 + e4) / (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8))

    cdef double den8 = sqrt(1 + e2 + e4 + e6 + e8)

    cdef double den9 = (1 - e2)**3 * (1 + 4 * e2 + e4)
    cdef double num9 = sqrt((1 - e2)**2 * (1 + 4 * e2 + e4) / (1 + 9 * e2 + 9 * e4 + e6))

    cdef double den10 = (1 - e2)**4 * (1 + e2) * (1 + e4)
    cdef double num10 = sqrt((1 - e2)**6 * (1 + e2) * (1 + e4) /
                             (1 + 4 * e2 + 10 * e4 + 20 * e6 + 10 * e8 +
                              4 * e10 + e12))

    cdef double den11 = sqrt(1 + e2 + e4 + e6 + e8 + e10)
    cdef double den12 = (1 - e2)**3

    # Parameters in loop
    cdef double x2, y2, x4, y4, xy, r2, r4, temp, x_c, y_c
    cdef int ii

    if (axis == "dx"):

        for ii in range(n):
            x_c = x[ii]
            y_c = y[ii]

            x2 = x_c * x_c
            y2 = y_c * y_c
            x4 = x2 * x2
            y4 = y2 * y2
            xy = x_c * y_c
            r2 = x2 + y2

            temp = Z[0] * 0 * x_c  # to make d an array with the same size as x

            temp += Z[1] * 2 * 1 / den1
            temp += Z[2] * 2 * 0

            temp += Z[3] * sqrt_3 * 4 * x_c / den2

            temp += Z[4] * sqrt_6 * 2 * y_c / den3
            temp += Z[5] * sqrt_6 * 2 * x_c / den3

            temp += Z[6] * sqrt_8 * 6 * xy * (1 + e2) / den4
            temp += Z[7] * sqrt_8 * ((9 * x2 + 3 * y2 - 2) * (1 + e2) - 2 * e4) / den4

            temp += Z[8] * sqrt_8 * 6 * xy / den5
            temp += Z[9] * sqrt_8 * (3 * x2 - 3 * y2) / den5

            temp += Z[10] * sqrt_5 * 12 * x_c * (2 * r2 - 1 - e2) / den6

            temp += Z[11] * sqrt_10 * (x_c * (16 * x2 - 6) *
                                       (1 + e2 + e4) - 6 * x_c * e6) * num7 / den7
            temp += Z[12] * sqrt_10 * (y_c * (24 * x2 + 8 * y2 - 6) *
                                       (1 + e2 + e4) - 6 * y_c * e6) * num7 / den7

            temp += Z[13] * sqrt_10 * 4 * x_c * (x2 - 3 * y2) / den8
            temp += Z[14] * sqrt_10 * 4 * y_c * (3 * x2 - y2) / den8

            temp += Z[15] * sqrt_12 * (
                3 * e8 - 36 * e6 * x2 - 12 * e6 * y2 + 12 * e6 +
                50 * e4 * x4 + 60 * e4 * x2 * y2 - 144 * e4 * x2 +
                10 * e4 * y4 - 48 * e4 * y2 + 30 * e4 + 200 * e2 * x4 + 240 *
                e2 * x2 * y2 - 144 * e2 * x2 + 40 * e2 * y4 - 48 * e2 * y2 +
                12 * e2 + 50 * x4 + 60 * x2 * y2 - 36 * x2 +
                10 * y4 - 12 * y2 + 3) * num9 / den9
            temp += Z[16] * sqrt_12 * (
                8 * xy * (5 * r2 * (1 + 4 * e2 + e4) -
                          (3 + 12 * e2 + 12 * e4 + 3 * e6))) * num9 / den9

            temp += Z[17] * sqrt_12 * (
                25 * (e6 + e4 + e2 + 1) * x4 +
                (- 12 * e8 - 30 * e6 * y2 - 12 * e6 - 30 * e4 * y2 - 12 * e4 -
                 30 * e2 * y2 - 12 * e2 - 30 * y2 - 12) * x2 + 12 * e8 * y2 -
                15 * e6 * y4 + 12 * e6 * y2 - 15 * e4 * y4 + 12 * e4 * y2 -
                15 * e2 * y4 + 12 * e2 * y2 - 15 * y4 + 12 * y2) * num10 / den10
            temp += Z[18] * sqrt_12 * (
                4.0 * xy * (15 * (e6 + e4 + e2 + 1) * x2 - 6 * e8 + 5 * e6 * y2 -
                            6 * e6 + 5 * e4 * y2 - 6 * e4 + 5 * e2 * y2 -
                            6 * e2 + 5 * y2 - 6)) * num10 / den10

            temp += Z[19] * sqrt_12 * 5 * (x2 * (x2 - 6 * y2) + y4) / den11
            temp += Z[20] * sqrt_12 * 20 * xy * (x2 - y2) / den11

            temp += Z[21] * sqrt_7 * 24 * x_c * (
                e4 - e2 * (5 * y2 - 3) + 5 * x4 - 5 * y2 + 5 * y4 -
                x2 * (5 * e2 - 10 * y2 + 5) + 1) / den12

            d[ii] = temp

    elif (axis == "dy"):

        for ii in range(n):

            x_c = x[ii]
            y_c = y[ii]

            x2 = x_c * x_c
            y2 = y_c * y_c
            x4 = x2 * x2
            y4 = y2 * y2
            xy = x_c * y_c
            r2 = x2 + y2

            temp = Z[0] * 0 * x_c   # to make d an array with the same size as x

            temp += Z[1] * 2 * 0
            temp += Z[2] * 2 * 1 / den1

            temp += Z[3] * sqrt_3 * 4 * y_c / den2

            temp += Z[4] * sqrt_6 * 2 * x_c / den3
            temp += Z[5] * sqrt_6 * (-2) * y_c / den3

            temp += Z[6] * sqrt_8 * ((1 + e2) *
                                     (3 * x2 + 9 * y2 - 2) - 2 * e4) / den4
            temp += Z[7] * sqrt_8 * 6 * xy * (1 + e2) / den4

            temp += Z[8] * sqrt_8 * (3 * x2 - 3 * y2) / den5
            temp += Z[9] * sqrt_8 * (-6) * xy / den5

            temp += Z[10] * sqrt_5 * 12 * y_c * (2 * r2 - 1 - e2) / den6

            temp += Z[11] * sqrt_10 * (y_c * (6 - 16 * y2) *
                                       (1 + e2 + e4) + 6 * y_c * e6) * num7 / den7
            temp += Z[12] * sqrt_10 * (x_c * (8 * x2 + 24 * y2 - 6) *
                                       (1 + e2 + e4) - 6 * x_c * e6) * num7 / den7

            temp += Z[13] * sqrt_10 * 4 * y_c * (y2 - 3 * x2) / den8
            temp += Z[14] * sqrt_10 * 4 * x_c * (x2 - 3 * y2) / den8

            temp += Z[15] * sqrt_12 * (
                -x_c * (24 * y_c + 4 * e2 * (24 * y_c - 40 * y_c * r2) +
                        2 * e4 * (48 * y_c - 20 * y_c * r2) + 24 * e6 * y_c -
                        40 * y_c * r2)) * num9 / den9
            temp += Z[16] * sqrt_12 * (
                3 * e8 - 12 * e6 * x2 - 36 * e6 * y2 + 12 * e6 + 10 * e4 * x4 +
                60 * e4 * x2 * y2 - 48 * e4 * x2 +
                50 * e4 * y4 - 144 * e4 * y2 + 30 * e4 + 40 * e2 * x4 + 240 *
                e2 * x2 * y2 - 48 * e2 * x2 + 200 * e2 * y4 - 144 * e2 * y2 +
                12 * e2 + 10 * x4 + 60 * x2 * y2 - 12 * x2 +
                50 * y4 - 36 * y2 + 3) * num9 / den9

            temp += Z[17] * sqrt_12 * (
                4.0 * xy * ((- 5) * (e6 + e4 + e2 + 1) * x2 + 6 * e8 -
                            15 * e6 * y2 + 6 * e6 - 15 * e4 * y2 +
                            6 * e4 - 15 * e2 * y2 + 6 * e2 -
                            15 * y2 + 6)) * num10 / den10
            temp += Z[18] * sqrt_12 * (
                - 12 * e8 * x2 + 12 * e8 * y2 + 15 * e6 * x4 +
                30 * e6 * x2 * y2 - 12 * e6 * x2 - 25 * e6 * y4 +
                12 * e6 * y2 + 15 * e4 * x4 + 30 * e4 * x2 * y2 - 12 * e4 * x2 -
                25 * e4 * y4 + 12 * e4 * y2 + 15 * e2 * x4 + 30 * e2 * x2 * y2 -
                12 * e2 * x2 - 25 * e2 * y4 + 12 * e2 * y2 + 15 * x4 +
                30 * x2 * y2 - 12 * x2 - 25 * y4 + 12 * y2) * num10 / den10

            temp += Z[19] * sqrt_12 * 20 * xy * (y2 - x2) / den11
            temp += Z[20] * sqrt_12 * 5 * (x2 * (x2 - 6 * y2) + y4) / den11

            temp += Z[21] * sqrt_7 * 24 * y_c * (
                e4 - e2 * (5 * x2 - 3) - 5 * x2 + 5 * x4 + 5 * y4 -
                y2 * (5 * e2 - 10 * x2 + 5) + 1) / den12

            d[ii] = temp

    elif (axis == "dx2"):

        for ii in range(n):

            x_c = x[ii]
            y_c = y[ii]

            x2 = x_c * x_c
            y2 = y_c * y_c
            x4 = x2 * x2
            y4 = y2 * y2
            xy = x_c * y_c
            r2 = x2 + y2
            r4 = r2 * r2

            temp = Z[0] * 0 * x_c  # to make d an array with the same size as x
            temp += Z[1] * 0
            temp += Z[2] * 0

            temp += Z[3] * sqrt_3 * 4 / den2
            temp += Z[4] * 0

            temp += Z[5] * sqrt_6 * 2 / den3

            temp += Z[6] * sqrt_8 * 6 * y_c * (1 + e2) / den4
            temp += Z[7] * sqrt_8 * 18 * x_c * (1 + e2) / den4

            temp += Z[8] * sqrt_8 * 6 * y_c / den5
            temp += Z[9] * sqrt_8 * 6 * x_c / den5

            temp += Z[10] * sqrt_5 * 12 * (6 * x2 + 2 * y2 - e2 - 1) / den6

            temp += Z[11] * sqrt_10 * ((48 * x2 - 6) *
                                       (1 + e2 + e4) - 6 * e6) * num7 / den7
            temp += Z[12] * sqrt_10 * 48 * xy * (1 + e2 + e4) * num7 / den7

            temp += Z[13] * sqrt_10 * 12 * (x2 - y2) / den8
            temp += Z[14] * sqrt_10 * 24 * xy / den8

            temp += Z[15] * sqrt_12 * (
                -8 * x_c * (9 * e6 - 25 * e4 * x2 - 15 * e4 * y2 + 36 * e4 -
                            100 * e2 * x2 - 60 * e2 * y2 + 36 * e2 - 25 * x2 -
                            15 * y2 + 9)) * num9 / den9
            temp += Z[16] * sqrt_12 * (
                -8 * y_c * (3 * e6 - 15 * e4 * x2 - 5 * e4 * y2 + 12 * e4 -
                            60 * e2 * x2 - 20 * e2 * y2 + 12 * e2 - 15 * x2 -
                            5 * y2 + 3)) * num9 / den9

            temp += Z[17] * sqrt_12 * (
                -4 * x_c * (6 * e8 - 25 * e6 * x2 + 15 * e6 * y2 + 6 * e6 -
                            25 * e4 * x2 + 15 * e4 * y2 + 6 * e4 - 25 * e2 * x2 +
                            15 * e2 * y2 + 6 * e2 - 25 * x2 +
                            15 * y2 + 6)) * num10 / den10
            temp += Z[18] * sqrt_12 * (
                -4 * y_c * (6 * e8 - 45 * e6 * x2 - 5 * e6 * y2 + 6 * e6 -
                            45 * e4 * x2 - 5 * e4 * y2 + 6 * e4 - 45 * e2 * x2 -
                            5 * e2 * y2 + 6 * e2 - 45 * x2 - 5 * y2 + 6)) * num10 / den10

            temp += Z[19] * sqrt_12 * 20 * x_c * (x2 - 3 * y2) / den11
            temp += Z[20] * sqrt_12 * 20 * y_c * (3 * x2 - y2) / den11

            temp += Z[21] * sqrt_7 * (
                480 * x2 * r2 + 120 * r4 + 24 * e4 - 360 * x2 - 120 * y2 -
                3 * e2 * (120 * x2 + 40 * y2 - 24) + 24) / den12

            d[ii] = temp

    elif (axis == "dy2"):

        for ii in range(n):

            x_c = x[ii]
            y_c = y[ii]

            x2 = x_c * x_c
            y2 = y_c * y_c
            x4 = x2 * x2
            y4 = y2 * y2
            xy = x_c * y_c
            r2 = x2 + y2
            r4 = r2 * r2

            temp = Z[0] * 0 * x_c  # to make d an array with the same size as x
            temp += Z[1] * 0
            temp += Z[2] * 0

            temp += Z[3] * sqrt_3 * 4 / den2
            temp += Z[4] * 0

            temp += Z[5] * sqrt_6 * (-2) / den3

            temp += Z[6] * sqrt_8 * (1 + e2) * 18 * y_c / den4
            temp += Z[7] * sqrt_8 * 6 * x_c * (1 + e2) / den4

            temp += Z[8] * sqrt_8 * (-6) * y_c / den5
            temp += Z[9] * sqrt_8 * (-6) * x_c / den5

            temp += Z[10] * sqrt_5 * 12 * (2 * x2 + 6 * y2 - e2 - 1) / den6

            temp += Z[11] * sqrt_10 * ((6 - 48 * y2) *
                                       (1 + e2 + e4) + 6 * e6) * num7 / den7
            temp += Z[12] * sqrt_10 * 48 * xy * (1 + e2 + e4) * num7 / den7

            temp += Z[13] * sqrt_10 * 12 * (y2 - x2) / den8
            temp += Z[14] * sqrt_10 * (-24) * xy / den8

            temp += Z[15] * sqrt_12 * (
                -8 * x_c * (3 * e6 - 5 * e4 * x2 - 15 * e4 * y2 + 12 * e4 -
                            20 * e2 * x2 - 60 * e2 * y2 + 12 * e2 - 5 * x2 -
                            15 * y2 + 3)) * num9 / den9
            temp += Z[16] * sqrt_12 * (
                -8 * y_c * (9 * e6 - 15 * e4 * x2 - 25 * e4 * y2 + 36 * e4 -
                            60 * e2 * x2 - 100 * e2 * y2 + 36 * e2 - 15 * x2 -
                            25 * y2 + 9)) * num9 / den9

            temp += Z[17] * sqrt_12 * (
                4 * x_c * (6 * e8 - 5 * e6 * x2 - 45 * e6 * y2 + 6 * e6 -
                           5 * e4 * x2 - 45 * e4 * y2 + 6 * e4 - 5 * e2 * x2 -
                           45 * e2 * y2 + 6 * e2 - 5 * x2 - 45 * y2 +
                           6)) * num10 / den10
            temp += Z[18] * sqrt_12 * (
                4 * y_c * (6 * e8 + 15 * e6 * x2 - 25 * e6 * y2 + 6 * e6 +
                           15 * e4 * x2 - 25 * e4 * y2 + 6 * e4 + 15 * e2 * x2 -
                           25 * e2 * y2 + 6 * e2 + 15 * x2 - 25 * y2 +
                           6)) * num10 / den10

            temp += Z[19] * sqrt_12 * 20 * x_c * (3 * y2 - x2) / den11
            temp += Z[20] * sqrt_12 * 20 * y_c * (y2 - 3 * x2) / den11

            temp += Z[21] * sqrt_7 * (
                480 * y2 * r2 + 120 * r4 + 24 * e4 - 120 * x2 - 360 * y2 -
                3 * e2 * (40 * x2 + 120 * y2 - 24) + 24) / den12

            d[ii] = temp

    elif (axis == "dxy"):

        for ii in range(n):

            x_c = x[ii]
            y_c = y[ii]

            x2 = x_c * x_c
            y2 = y_c * y_c
            x4 = x2 * x2
            y4 = y2 * y2
            xy = x_c * y_c
            r2 = x2 + y2
            r4 = r2 * r2

            temp = Z[0] * 0 * x_c  # to make d an array with the same size as x
            temp += Z[1] * 0
            temp += Z[2] * 0
            temp += Z[3] * 0

            temp += Z[4] * sqrt_6 * 2 / den3
            temp += Z[5] * 0

            temp += Z[6] * sqrt_8 * (1 + e2) * (6 * x_c) / den4
            temp += Z[7] * sqrt_8 * 6 * y_c * (1 + e2) / den4

            temp += Z[8] * sqrt_8 * 6 * x_c / den5
            temp += Z[9] * sqrt_8 * (-6) * y_c / den5

            temp += Z[10] * sqrt_5 * 48 * xy / den6

            temp += Z[11] * sqrt_10 * 0
            temp += Z[12] * sqrt_10 * ((24 * x2 + 24 * y2 - 6) *
                                       (1 + e2 + e4) - 6 * e6) * num7 / den7

            temp += Z[13] * sqrt_10 * (-24) * xy / den8
            temp += Z[14] * sqrt_10 * 12 * (x2 - y2) / den8

            temp += Z[15] * sqrt_12 * (
                -8 * y_c * (3 * e6 - 15 * e4 * x2 - 5 * e4 * y2 + 12 * e4 -
                            60 * e2 * x2 - 20 * e2 * y2 + 12 * e2 - 15 * x2 -
                            5 * y2 + 3)) * num9 / den9
            temp += Z[16] * sqrt_12 * (
                -8 * x_c * (3 * e6 - 5 * e4 * x2 - 15 * e4 * y2 + 12 * e4 -
                            20 * e2 * x2 - 60 * e2 * y2 + 12 * e2 - 5 * x2 -
                            15 * y2 + 3)) * num9 / den9

            temp += Z[17] * sqrt_12 * (
                12 * y_c * (2 * e8 - 5 * e6 * r2 + 2 * e6 - 5 * e4 * r2 + 2 * e4 -
                            5 * e2 * r2 + 2 * e2 - 5 * r2 + 2)) * num10 / den10
            temp += Z[18] * sqrt_12 * (
                -12 * x_c * (2 * e8 - 5 * e6 * r2 + 2 * e6 -
                             5 * e4 * r2 + 2 * e4 - 5 * e2 * r2 + 2 * e2 -
                             5 * r2 + 2)) * num10 / den10

            temp += Z[19] * sqrt_12 * 20 * y_c * (y2 - 3 * x2) / den11
            temp += Z[20] * sqrt_12 * 20 * x_c * (x2 - 3 * y2) / den11

            temp += Z[21] * sqrt_7 * 240 * xy * (2 * r2 - 1 - e2) / den12

            d[ii] = temp


def ZernikeAnnularEval(ndarray[np.float64_t, ndim=1] Z not None, ndarray[np.float64_t, ndim=1] x not None, 
                       ndarray[np.float64_t, ndim=1] y not None, e):

    cdef Py_ssize_t n = x.shape[0]
    cdef ndarray[np.float64_t, ndim=1] S = np.zeros_like(x)

    _ZernikeAnnularEval(&S[0], &Z[0], &x[0], &y[0], e, n)

    return S

cdef _ZernikeAnnularEval(double *S, double *Z, double *x, double *y, double e, int n):

    # Parameters of constant
    cdef double e2 = e**2
    cdef double e4 = e2*e2
    cdef double e6 = e4*e2
    cdef double e8 = e6*e2
    cdef double e10 = e8*e2
    cdef double e12 = e10*e2
    cdef double e14 = e12*e2

    cdef double sqrt_3 = sqrt(3)
    cdef double sqrt_5 = sqrt(5)
    cdef double sqrt_6 = sqrt(6)
    cdef double sqrt_7 = sqrt(7)
    cdef double sqrt_8 = sqrt(8)
    cdef double sqrt_10 = sqrt(10)
    cdef double sqrt_12 = sqrt(12)
    cdef double sqrt_14 = sqrt(14)

    cdef double den1 = sqrt(1 + e2)
    cdef double den2 = 1 - e2
    cdef double den3 = sqrt(1 + e2 + e4)
    cdef double den4 = sqrt((1 - e2)**2 * (1 + e2) * (1 + 4 * e2 + e4))
    cdef double den5 = sqrt(1 + e2 + e4 + e6)
    cdef double den6 = (1 - e2)**2

    cdef double den7 = (1 - e2)**3 * (1 + e2 + e4)
    cdef double num7 = sqrt((1 - e2)**4 * (1 + e2 + e4) /
                            (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8))

    cdef double den8 = sqrt(1 + e2 + e4 + e6 + e8)

    cdef double den9 = (1 - e2)**3 * (1 + 4 * e2 + e4)
    cdef double num9E = sqrt((1 - e2)**2 * (1 + 4 * e2 + e4) /
                             (1 + 9 * e2 + 9 * e4 + e6))

    cdef double den10 = (1 - e2)**4 * (1 + e2) * (1 + e4)
    cdef double num10E = sqrt((1 - e2)**6 * (1 + e2) * (1 + e4) /
                              (1 + 4 * e2 + 10 * e4 + 20 * e6 + 10 * e8 + 4 * e10 + e12))

    cdef double den11 = sqrt(1 + e2 + e4 + e6 + e8 + e10)
    cdef double den12 = (1 - e2)**3

    cdef double num11a = 15 * (1 + 4*e2 + 10*e4 + 4*e6 + e8)
    cdef double num11b = -20 * (1 + 4*e2 + 10*e4 + 10*e6 + 4*e8 + e10)
    cdef double num11c = 6 * (1 + 4*e2 + 10*e4 + 20*e6 + 10*e8 + 4*e10 + e12)
    cdef double den13 = (1-e2)**2 * sqrt((1 + 4*e2 + 10*e4 + 4*e6 + e8) * (
        1 + 9*e2 + 45*e4 + 65*e6 + 45*e8 + 9*e10 + e12))

    cdef double num12 = -5 * (1 - e12) / (1 - e10)
    cdef double den14 = sqrt(1 / (1-e2) * (36*(1-e14) - (35 * (1 - e12)**2) / (1 - e10)))

    cdef double num13 = sqrt((1 - e2) / (1 - e14))

    # Parameter in loop
    cdef double r, r2, r3, r4, r5, r6
    cdef double t, t2, t3, t4, t5, t6
    cdef double s, s2, s3, s4, s5, s6
    cdef double c, c2, c3, c4, c5, c6

    cdef double temp, numQ, x_c, y_c, Rnl
    cdef int ii

    for ii in range(n):

        x_c = x[ii]
        y_c = y[ii]

        r2 = x_c**2 + y_c**2
        r = sqrt(r2)
        r3 = r2 * r
        r4 = r2 * r2
        r5 = r3 * r2
        r6 = r3 * r3

        t = atan2(y_c, x_c)
        s = sin(t)
        c = cos(t)

        t2 = 2*t
        t3 = 3*t
        t4 = 4*t
        t5 = 5*t
        t6 = 6*t

        s2 = sin(t2)
        c2 = cos(t2)
        s3 = sin(t3)
        c3 = cos(t3)
        s4 = sin(t4)
        c4 = cos(t4)
        s5 = sin(t5)
        c5 = cos(t5)
        s6 = sin(t6)
        c6 = cos(t6)

        temp = Z[0] * (1 + 0 * x_c)

        Rnl = 2 * r / den1
        temp += Z[1] * Rnl * c
        temp += Z[2] * Rnl * s

        temp += Z[3] * sqrt_3 * (2 * r2 - 1 - e2) / den2

        Rnl = sqrt_6 * r2 / den3
        temp += Z[4] * Rnl * s2
        temp += Z[5] * Rnl * c2

        Rnl = sqrt_8 * (3 * r3 - 2 * r - 2 * e4 * r + e2 * r * (3 * r2 - 2)) / den4
        temp += Z[6] * Rnl * s
        temp += Z[7] * Rnl * c

        Rnl = sqrt_8 * r3 / den5
        temp += Z[8] * Rnl * s3
        temp += Z[9] * Rnl * c3

        temp += Z[10] * sqrt_5 * (6 * r4 - 6 * r2 + 1 +
                                  e4 + e2 * (4 - 6 * r2)) / den6

        Rnl = sqrt_10 * (4 * r4 - 3 * r2 - 3 * e6 * r2 - e2 * r2 * (3 - 4 * r2) -
                         e4 * r2 * (3 - 4 * r2)) * num7 / den7
        temp += Z[11] * Rnl * c2
        temp += Z[12] * Rnl * s2

        Rnl = sqrt_10 * r4 / den8
        temp += Z[13] * Rnl * c4
        temp += Z[14] * Rnl * s4

        numQ = 10 * r5 - 12 * r3 + 3 * r + 3 * e8 * r - 12 * e6 * r * (r2 - 1) + \
            2 * e4 * r * (15 - 24 * r2 + 5 * r4) + 4 * e2 * r * (3 - 12 * r2 + 10 * r4)
        Rnl = sqrt_12 * num9E * numQ / den9
        temp += Z[15] * Rnl * c
        temp += Z[16] * Rnl * s

        numQ = r3 * (5 * r2 - 4 - 4 * e8 - e2 * (4 - 5 * r2) -
                     e4 * (4 - 5 * r2) - e6 * (4 - 5 * r2))
        Rnl = sqrt_12 * num10E * numQ / den10
        temp += Z[17] * Rnl * c3
        temp += Z[18] * Rnl * s3

        Rnl = sqrt_12 * r5 / den11
        temp += Z[19] * Rnl * c5
        temp += Z[20] * Rnl * s5

        temp += Z[21] * sqrt_7 * (
            20 * r6 - 30 * r4 + 12 * r2 - 1 - e6 +
            3 * e4 * (-3 + 4 * r2) - 3 * e2 * (3 - 12 * r2 + 10 * r4)) / den12

        Rnl = sqrt_14 * (num11a*r6 + num11b*r4 + num11c*r2) / den13
        temp += Z[22] * Rnl * s2
        temp += Z[23] * Rnl * c2

        Rnl = sqrt_14 * (6*r6 + num12*r4) / den14
        temp += Z[24] * Rnl * s4
        temp += Z[25] * Rnl * c4

        Rnl = sqrt_14 * num13 * r6
        temp += Z[26] * Rnl * s6
        temp += Z[27] * Rnl * c6

        S[ii] = temp


def poly10_2D(ndarray[np.float64_t, ndim=1] c not None, ndarray[np.float64_t,
              ndim=1] x not None, ndarray[np.float64_t, ndim=1] y not None):

    cdef Py_ssize_t n = x.shape[0]
    cdef ndarray[np.float64_t, ndim=1] cy_out = np.zeros_like(x)

    _poly10_2D(&cy_out[0], &c[0], &x[0], &y[0], n)

    return cy_out


cdef _poly10_2D(double *cy_out, double *c, double *x, double *y, int n):

    cdef double x_c, y_c
    cdef int ii

    for ii in range(n):
        x_c = x[ii]
        y_c = y[ii]
        cy_out[ii] = c[0] + c[1] * x_c + c[2] * y_c + c[3] * x_c * x_c + \
            c[4] * x_c * y_c + c[5] * y_c * y_c + c[6] * x_c**3 + \
            c[7] * x_c**2 * y_c + c[8] * x_c * y_c**2 + c[9] * y_c**3 + \
            c[10] * x_c**4 + c[11] * x_c**3 * y_c + c[12] * x_c**2 * y_c**2 + \
            c[13] * x_c * y_c**3 + c[14] * y_c**4 + c[15] * x_c**5 + \
            c[16] * x_c**4 * y_c + c[17] * x_c**3 * y_c**2 + c[18] * x_c**2 * y_c**3 + \
            c[19] * x_c * y_c**4 + c[20] * y_c**5 + c[21] * x_c**6 + \
            c[22] * x_c**5 * y_c + c[23] * x_c**4 * y_c**2 + c[24] * x_c**3 * y_c**3 + \
            c[25] * x_c**2 * y_c**4 + c[26] * x_c * y_c**5 + c[27] * y_c**6 + \
            c[28] * x_c**7 + c[29] * x_c**6 * y_c + c[30] * x_c**5 * y_c**2 + \
            c[31] * x_c**4 * y_c**3 + c[32] * x_c**3 * y_c**4 + c[33] * x_c**2 * y_c**5 + \
            c[34] * x_c * y_c**6 + c[35] * y_c**7 + c[36] * x_c**8 + \
            c[37] * x_c**7 * y_c + c[38] * x_c**6 * y_c**2 + c[39] * x_c**5 * y_c**3 + \
            c[40] * x_c**4 * y_c**4 + c[41] * x_c**3 * y_c**5 + c[42] * x_c**2 * y_c**6 + \
            c[43] * x_c * y_c**7 + c[44] * y_c**8 + c[45] * x_c**9 + \
            c[46] * x_c**8 * y_c + c[47] * x_c**7 * y_c**2 + c[48] * x_c**6 * y_c**3 + \
            c[49] * x_c**5 * y_c**4 + c[50] * x_c**4 * y_c**5 + c[51] * x_c**3 * y_c**6 + \
            c[52] * x_c**2 * y_c**7 + c[53] * x_c * y_c**8 + c[54] * y_c**9 + \
            c[55] * x_c**10 + c[56] * x_c**9 * y_c + c[57] * x_c**8 * y_c**2 + \
            c[58] * x_c**7 * y_c**3 + c[59] * x_c**6 * y_c**4 + c[60] * x_c**5 * y_c**5 + \
            c[61] * x_c**4 * y_c**6 + c[62] * x_c**3 * y_c**7 + c[63] * x_c**2 * y_c**8 + \
            c[64] * x_c * y_c**9 + c[65] * y_c**10


def poly10Grad(ndarray[np.float64_t, ndim=1] c not None, ndarray[np.float64_t,
               ndim=1] x not None, ndarray[np.float64_t, ndim=1] y not None, axis):

    cdef Py_ssize_t n = x.shape[0]
    cdef ndarray[np.float64_t, ndim=1] cy_out = np.zeros_like(x)

    _poly10Grad(&cy_out[0], &c[0], &x[0], &y[0], axis, n)

    return cy_out


cdef _poly10Grad(double *cy_out, double *c, double *x, double *y, str axis, int n):

    cdef double x_c, y_c
    cdef int ii

    if (axis == "dx"):
        for ii in range(n):
            x_c = x[ii]
            y_c = y[ii]
            cy_out[ii] = c[1] + c[3] * 2 * x_c + c[4] * y_c + c[6] * 3 * x_c**2 + \
                c[7] * 2 * x_c * y_c + c[8] * y_c**2 + c[10] * 4 * x_c**3 + \
                c[11] * 3 * x_c**2 * y_c + c[12] * 2 * x_c * y_c**2 + c[13] * y_c**3 + \
                c[15] * 5 * x_c**4 + c[16] * 4 * x_c**3 * y_c + \
                c[17] * 3 * x_c**2 * y_c**2 + \
                c[18] * 2 * x_c * y_c**3 + c[19] * y_c**4 + c[21] * 6 * x_c**5 + \
                c[22] * 5 * x_c**4 * y_c + c[23] * 4 * x_c**3 * y_c**2 + \
                c[24] * 3 * x_c**2 * y_c**3 + c[25] * 2 * x_c * y_c**4 + c[26] * y_c**5 + \
                c[28] * 7 * x_c**6 + c[29] * 6 * x_c**5 * y_c + \
                c[30] * 5 * x_c**4 * y_c**2 + \
                c[31] * 4 * x_c**3 * y_c**3 + c[32] * 3 * x_c**2 * y_c**4 + \
                c[33] * 2 * x_c * y_c**5 + c[34] * y_c**6 + \
                c[36] * 8 * x_c**7 + c[37] * 7 * x_c**6 * y_c + \
                c[38] * 6 * x_c**5 * y_c**2 + \
                c[39] * 5 * x_c**4 * y_c**3 + c[40] * 4 * x_c**3 * y_c**4 + \
                c[41] * 3 * x_c**2 * y_c**5 + c[42] * 2 * x_c * y_c**6 + c[43] * y_c**7 + \
                c[45] * 9 * x_c**8 + c[46] * 8 * x_c**7 * y_c + \
                c[47] * 7 * x_c**6 * y_c**2 + \
                c[48] * 6 * x_c**5 * y_c**3 + c[49] * 5 * x_c**4 * y_c**4 + \
                c[50] * 4 * x_c**3 * y_c**5 + c[51] * 3 * x_c**2 * y_c**6 + \
                c[52] * 2 * x_c * y_c**7 + c[53] * y_c**8 + c[55] * 10 * x_c**9 + \
                c[56] * 9 * x_c**8 * y_c + c[57] * 8 * x_c**7 * y_c**2 + \
                c[58] * 7 * x_c**6 * y_c**3 + c[59] * 6 * x_c**5 * y_c**4 + \
                c[60] * 5 * x_c**4 * y_c**5 + c[61] * 4 * x_c**3 * y_c**6 + \
                c[62] * 3 * x_c**2 * y_c**7 + c[63] * 2 * x_c * y_c**8 + c[64] * y_c**9
    elif (axis == "dy"):
        for ii in range(n):
            x_c = x[ii]
            y_c = y[ii]
            cy_out[ii] = c[2] + c[4] * x_c + c[5] * 2 * y_c + c[7] * x_c**2 + \
                c[8] * x_c * 2 * y_c + c[9] * 3 * y_c**2 + c[11] * x_c**3 + \
                c[12] * x_c**2 * 2 * y_c + c[13] * x_c * 3 * y_c**2 + c[14] * 4 * y_c**3 + \
                c[16] * x_c**4 + c[17] * x_c**3 * 2 * y_c + c[18] * x_c**2 * 3 * y_c**2 + \
                c[19] * x_c * 4 * y_c**3 + c[20] * 5 * y_c**4 + c[22] * x_c**5 + \
                c[23] * x_c**4 * 2 * y_c + c[24] * x_c**3 * 3 * y_c**2 + \
                c[25] * x_c**2 * 4 * y_c**3 + c[26] * x_c * 5 * y_c**4 + \
                c[27] * 6 * y_c**5 + c[29] * x_c**6 + c[30] * x_c**5 * 2 * y_c + \
                c[31] * x_c**4 * 3 * y_c**2 + c[32] * x_c**3 * 4 * y_c**3 + \
                c[33] * x_c**2 * 5 * y_c**4 + c[34] * x_c * 6 * y_c**5 + \
                c[35] * 7 * y_c**6 + \
                c[37] * x_c**7 + c[38] * x_c**6 * 2 * y_c + c[39] * x_c**5 * 3 * y_c**2 + \
                c[40] * x_c**4 * 4 * y_c**3 + c[41] * x_c**3 * 5 * y_c**4 + \
                c[42] * x_c**2 * 6 * y_c**5 + c[43] * x_c * 7 * y_c**6 + \
                c[44] * 8 * y_c**7 + \
                c[46] * x_c**8 + c[47] * x_c**7 * 2 * y_c + c[48] * x_c**6 * 3 * y_c**2 + \
                c[49] * x_c**5 * 4 * y_c**3 + c[50] * x_c**4 * 5 * y_c**4 + \
                c[51] * x_c**3 * 6 * y_c**5 + c[52] * x_c**2 * 7 * y_c**6 + \
                c[53] * x_c * 8 * y_c**7 + c[54] * 9 * y_c**8 + c[56] * x_c**9 + \
                c[57] * x_c**8 * 2 * y_c + c[58] * x_c**7 * 3 * y_c**2 + \
                c[59] * x_c**6 * 4 * y_c**3 + c[60] * x_c**5 * 5 * y_c**4 + \
                c[61] * x_c**4 * 6 * y_c**5 + c[62] * x_c**3 * 7 * y_c**6 + \
                c[63] * x_c**2 * 8 * y_c**7 + c[64] * x_c * 9 * y_c**8 + c[65] * 10 * y_c**9
