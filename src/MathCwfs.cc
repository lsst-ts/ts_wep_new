/*
 * This file is part of ts_wep.
 *
 * Developed for the LSST Telescope and Site Systems.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cmath>

#include <pybind11/numpy.h>

#include "lsst/ts/wep/cwfs/MathCwfs.h"

namespace py = pybind11;

namespace lsst {
namespace ts {
namespace wep {
namespace cwfs {

struct ArrayNpInfo {
    size_t size;
    double *buf;
};

ArrayNpInfo getNpArrayInfo(py::array arrayNp) {

    ArrayNpInfo info;

    py::buffer_info bufInfo = arrayNp.request();
    if (bufInfo.ndim != 1) {
        throw std::runtime_error("Number of dimensions must be one.");
    }

    info.size = bufInfo.size;
    info.buf = (double *)bufInfo.ptr;

    return info;
}

py::array_t<double> zernikeAnnularEval(py::array_t<double> arrayZk,
                                       py::array_t<double> arrayX,
                                       py::array_t<double> arrayY, double e) {
    // Get the numpy array information
    ArrayNpInfo infoZk = getNpArrayInfo(arrayZk);
    ArrayNpInfo infoX = getNpArrayInfo(arrayX);
    ArrayNpInfo infoY = getNpArrayInfo(arrayY);

    // No pointer is passed, so NumPy will allocate the buffer
    size_t n = infoX.size;
    auto result = py::array_t<double>(n);
    ArrayNpInfo infoResult = getNpArrayInfo(result);

    // Assign the variables
    double *Z = infoZk.buf;
    double *x = infoX.buf;
    double *y = infoY.buf;
    double *S = infoResult.buf;

    // Parameters of constant
    double e2 = pow(e, 2);
    double e4 = e2 * e2;
    double e6 = e4 * e2;
    double e8 = e6 * e2;
    double e10 = e8 * e2;
    double e12 = e10 * e2;
    double e14 = e12 * e2;

    double sqrt_3 = sqrt(3);
    double sqrt_5 = sqrt(5);
    double sqrt_6 = sqrt(6);
    double sqrt_7 = sqrt(7);
    double sqrt_8 = sqrt(8);
    double sqrt_10 = sqrt(10);
    double sqrt_12 = sqrt(12);
    double sqrt_14 = sqrt(14);

    double den1 = sqrt(1 + e2);
    double den2 = 1 - e2;
    double den3 = sqrt(1 + e2 + e4);
    double den4 = sqrt(pow(1 - e2, 2) * (1 + e2) * (1 + 4 * e2 + e4));
    double den5 = sqrt(1 + e2 + e4 + e6);
    double den6 = pow(1 - e2, 2);

    double den7 = pow(1 - e2, 3) * (1 + e2 + e4);
    double num7 = sqrt(pow(1 - e2, 4) * (1 + e2 + e4) /
                       (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8));

    double den8 = sqrt(1 + e2 + e4 + e6 + e8);

    double den9 = pow(1 - e2, 3) * (1 + 4 * e2 + e4);
    double num9E =
        sqrt(pow(1 - e2, 2) * (1 + 4 * e2 + e4) / (1 + 9 * e2 + 9 * e4 + e6));

    double den10 = pow(1 - e2, 4) * (1 + e2) * (1 + e4);
    double num10E =
        sqrt(pow(1 - e2, 6) * (1 + e2) * (1 + e4) /
             (1 + 4 * e2 + 10 * e4 + 20 * e6 + 10 * e8 + 4 * e10 + e12));

    double den11 = sqrt(1 + e2 + e4 + e6 + e8 + e10);
    double den12 = pow(1 - e2, 3);

    double num11a = 15 * (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8);
    double num11b = -20 * (1 + 4 * e2 + 10 * e4 + 10 * e6 + 4 * e8 + e10);
    double num11c =
        6 * (1 + 4 * e2 + 10 * e4 + 20 * e6 + 10 * e8 + 4 * e10 + e12);
    double den13 =
        pow(1 - e2, 2) *
        sqrt((1 + 4 * e2 + 10 * e4 + 4 * e6 + e8) *
             (1 + 9 * e2 + 45 * e4 + 65 * e6 + 45 * e8 + 9 * e10 + e12));

    double num12 = -5 * (1 - e12) / (1 - e10);
    double den14 = sqrt(1 / (1 - e2) *
                        (36 * (1 - e14) - (35 * pow(1 - e12, 2)) / (1 - e10)));

    double num13 = sqrt((1 - e2) / (1 - e14));

    // Parameter in loop
    double r, r2, r3, r4, r5, r6;
    double t, t2, t3, t4, t5, t6;
    double s, s2, s3, s4, s5, s6;
    double c, c2, c3, c4, c5, c6;

    double temp, numQ, x_c, y_c, Rnl;

    for (size_t ii = 0; ii < n; ii++) {
        x_c = x[ii];
        y_c = y[ii];

        r2 = pow(x_c, 2) + pow(y_c, 2);
        r = sqrt(r2);
        r3 = r2 * r;
        r4 = r2 * r2;
        r5 = r3 * r2;
        r6 = r3 * r3;

        t = atan2(y_c, x_c);
        s = sin(t);
        c = cos(t);

        t2 = 2 * t;
        t3 = 3 * t;
        t4 = 4 * t;
        t5 = 5 * t;
        t6 = 6 * t;

        s2 = sin(t2);
        c2 = cos(t2);
        s3 = sin(t3);
        c3 = cos(t3);
        s4 = sin(t4);
        c4 = cos(t4);
        s5 = sin(t5);
        c5 = cos(t5);
        s6 = sin(t6);
        c6 = cos(t6);

        temp = Z[0] * (1 + 0 * x_c);

        Rnl = 2 * r / den1;
        temp += Z[1] * Rnl * c;
        temp += Z[2] * Rnl * s;

        temp += Z[3] * sqrt_3 * (2 * r2 - 1 - e2) / den2;

        Rnl = sqrt_6 * r2 / den3;
        temp += Z[4] * Rnl * s2;
        temp += Z[5] * Rnl * c2;

        Rnl = sqrt_8 * (3 * r3 - 2 * r - 2 * e4 * r + e2 * r * (3 * r2 - 2)) /
              den4;
        temp += Z[6] * Rnl * s;
        temp += Z[7] * Rnl * c;

        Rnl = sqrt_8 * r3 / den5;
        temp += Z[8] * Rnl * s3;
        temp += Z[9] * Rnl * c3;

        temp += Z[10] * sqrt_5 *
                (6 * r4 - 6 * r2 + 1 + e4 + e2 * (4 - 6 * r2)) / den6;

        Rnl = sqrt_10 *
              (4 * r4 - 3 * r2 - 3 * e6 * r2 - e2 * r2 * (3 - 4 * r2) -
               e4 * r2 * (3 - 4 * r2)) *
              num7 / den7;
        temp += Z[11] * Rnl * c2;
        temp += Z[12] * Rnl * s2;

        Rnl = sqrt_10 * r4 / den8;
        temp += Z[13] * Rnl * c4;
        temp += Z[14] * Rnl * s4;

        numQ = 10 * r5 - 12 * r3 + 3 * r + 3 * e8 * r - 12 * e6 * r * (r2 - 1) +
               2 * e4 * r * (15 - 24 * r2 + 5 * r4) +
               4 * e2 * r * (3 - 12 * r2 + 10 * r4);
        Rnl = sqrt_12 * num9E * numQ / den9;
        temp += Z[15] * Rnl * c;
        temp += Z[16] * Rnl * s;

        numQ = r3 * (5 * r2 - 4 - 4 * e8 - e2 * (4 - 5 * r2) -
                     e4 * (4 - 5 * r2) - e6 * (4 - 5 * r2));
        Rnl = sqrt_12 * num10E * numQ / den10;
        temp += Z[17] * Rnl * c3;
        temp += Z[18] * Rnl * s3;

        Rnl = sqrt_12 * r5 / den11;
        temp += Z[19] * Rnl * c5;
        temp += Z[20] * Rnl * s5;

        temp += Z[21] * sqrt_7 *
                (20 * r6 - 30 * r4 + 12 * r2 - 1 - e6 + 3 * e4 * (-3 + 4 * r2) -
                 3 * e2 * (3 - 12 * r2 + 10 * r4)) /
                den12;

        Rnl = sqrt_14 * (num11a * r6 + num11b * r4 + num11c * r2) / den13;
        temp += Z[22] * Rnl * s2;
        temp += Z[23] * Rnl * c2;

        Rnl = sqrt_14 * (6 * r6 + num12 * r4) / den14;
        temp += Z[24] * Rnl * s4;
        temp += Z[25] * Rnl * c4;

        Rnl = sqrt_14 * num13 * r6;
        temp += Z[26] * Rnl * s6;
        temp += Z[27] * Rnl * c6;

        S[ii] = temp;
    }

    return result;
}

py::array_t<double> zernikeAnnularJacobian(py::array_t<double> arrayZk,
                                           py::array_t<double> arrayX,
                                           py::array_t<double> arrayY, double e,
                                           std::string atype) {
    // Get the numpy array information
    ArrayNpInfo infoZk = getNpArrayInfo(arrayZk);
    ArrayNpInfo infoX = getNpArrayInfo(arrayX);
    ArrayNpInfo infoY = getNpArrayInfo(arrayY);

    // No pointer is passed, so NumPy will allocate the buffer
    size_t n = infoX.size;
    auto result = py::array_t<double>(n);
    ArrayNpInfo infoResult = getNpArrayInfo(result);

    // Assign the variables
    double *Z = infoZk.buf;
    double *x = infoX.buf;
    double *y = infoY.buf;
    double *out = infoResult.buf;

    // Parameters of constant
    double e2 = pow(e, 2);
    double e4 = e2 * e2;
    double e6 = e4 * e2;
    double e8 = e6 * e2;
    double e10 = e8 * e2;
    double e12 = e10 * e2;
    double e14 = e12 * e2;
    double e16 = e14 * e2;

    double sqrt_3 = sqrt(3);
    double sqrt_5 = sqrt(5);
    double sqrt_6 = sqrt(6);
    double sqrt_7 = sqrt(7);
    double sqrt_8 = sqrt(8);
    double sqrt_10 = sqrt(10);
    double sqrt_12 = sqrt(12);

    // 1st order
    double den1 = 1 - e2;
    double den2 = sqrt(pow(1 - e2, 2) * (1 + e2) * (1 + 4 * e2 + e4));
    double den3 = pow(1 - e2, 2);

    double den4 = pow(1 - e2, 3) * (1 + e2 + e4);
    double num4 = sqrt(pow(1 - e2, 4) * (1 + e2 + e4) /
                       (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8));

    double den5 = pow(1 - e2, 3) * (1 + 4 * e2 + e4);
    double num5 =
        sqrt(pow(1 - e2, 2) * (1 + 4 * e2 + e4) / (1 + 9 * e2 + 9 * e4 + e6));

    double den6 = pow(1 - e2, 4) * (1 + e2) * (1 + e4);
    double num6 =
        sqrt(pow(1 - e2, 6) * (1 + e2) * (1 + e4) /
             (1 + 4 * e2 + 10 * e4 + 20 * e6 + 10 * e8 + 4 * e10 + e12));

    double den7 = pow(1 - e2, 3);

    // 2nd order
    double den2_2 = (1 + e2 + e4);
    double den2_3 = pow(1 - e2, 2) * (1 + e2) * (1 + 4 * e2 + e4);
    double den2_4 = (1 + e2 + e4 + e6);
    double den2_5 = pow(1 - e2, 4);

    double den2_6 = pow(1 - e2, 6) * pow(1 + e2 + e4, 2);
    double num2_6 =
        (pow(1 - e2, 4) * (1 + e2 + e4) / (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8));

    double den2_7 = (1 + e2 + e4 + e6 + e8);

    double den2_8 = pow(1 - e2, 6) * pow(1 + 4 * e2 + e4, 2);
    double num2_8 =
        pow(1 - e2, 2) * (1 + 4 * e2 + e4) / (1 + 9 * e2 + 9 * e4 + e6);

    double den2_9 = pow(1 - e2, 8) * pow(1 + e2, 2) * pow(1 + e4, 2);
    double num2_9 = pow(1 - e2, 6) * (1 + e2) * (1 + e4) /
                    (1 + 4 * e2 + 10 * e4 + 20 * e6 + 10 * e8 + 4 * e10 + e12);

    double den2_10 = (1 + e2 + e4 + e6 + e8 + e10);
    double den2_11 = pow(1 - e2, 6);

    // Parameters in loop
    double x2, y2, x4, y4, xy, r2, x6, y6, temp, x_c, y_c;

    if (atype == "1st") {
        for (size_t ii = 0; ii < n; ii++) {
            x_c = x[ii];
            y_c = y[ii];

            x2 = x_c * x_c;
            y2 = y_c * y_c;
            xy = x_c * y_c;
            r2 = x2 + y2;
            x4 = x2 * x2;
            y4 = y2 * y2;

            // to make d an array with the same size as x
            temp = Z[0] * 0 * x_c;
            temp += Z[1] * 0;
            temp += Z[2] * 0;

            temp += Z[3] * sqrt_3 * 8 / den1;
            temp += Z[4] * sqrt_6 * 0;
            temp += Z[5] * sqrt_6 * 0;

            temp += Z[6] * sqrt_8 * 24 * y_c * (1 + e2) / den2;
            temp += Z[7] * sqrt_8 * 24 * x_c * (1 + e2) / den2;
            temp += Z[8] * sqrt_8 * 0;
            temp += Z[9] * sqrt_8 * 0;

            temp += Z[10] * sqrt_5 * (96 * r2 - 24 * (1 + e2)) / den3;

            temp +=
                Z[11] * sqrt_10 * 48 * (x2 - y2) * (1 + e2 + e4) * num4 / den4;
            temp += Z[12] * sqrt_10 * 96 * xy * (1 + e2 + e4) * num4 / den4;
            temp += Z[13] * sqrt_10 * 0;
            temp += Z[14] * sqrt_10 * 0;

            temp +=
                Z[15] * sqrt_12 * 48 * x_c *
                (5 * r2 * (1 + 4 * e2 + e4) - 2 * (1 + 4 * e2 + 4 * e4 + e6)) *
                num5 / den5;
            temp +=
                Z[16] * sqrt_12 * 48 * y_c *
                (5 * r2 * (1 + 4 * e2 + e4) - 2 * (1 + 4 * e2 + 4 * e4 + e6)) *
                num5 / den5;

            temp += Z[17] * sqrt_12 * 80.0 * x_c * (x2 - 3.0 * y2) * (1 + e2) *
                    (1 + e4) * num6 / den6;
            temp += Z[18] * sqrt_12 * 80.0 * y_c * (3 * x2 - y2) * (1 + e2) *
                    (1 + e4) * num6 / den6;
            temp += Z[19] * sqrt_12 * 0;
            temp += Z[20] * sqrt_12 * 0;

            temp += Z[21] * sqrt_7 * 48 *
                    (e4 - 10 * e2 * x2 - 10 * e2 * y2 + 3 * e2 + 15 * x4 +
                     30 * x2 * y2 - 10 * x2 + 15 * y4 - 10 * y2 + 1) /
                    den7;

            out[ii] = temp;
        }
    } else if (atype == "2nd") {
        for (size_t ii = 0; ii < n; ii++) {
            x_c = x[ii];
            y_c = y[ii];

            x2 = x_c * x_c;
            y2 = y_c * y_c;
            xy = x_c * y_c;
            r2 = x2 + y2;
            x4 = x2 * x2;
            x6 = x4 * x2;
            y4 = y2 * y2;
            y6 = y4 * y2;

            // to make d an array with the same size as x
            temp = pow(Z[0], 2) * 0 * x_c;
            temp += pow(Z[1], 2) * 0;
            temp += pow(Z[2], 2) * 0;

            temp += pow(Z[3], 2) * (3) * 16 / den1 / den1;

            temp += pow(Z[4], 2) * (6) * (-4) / den2_2;
            temp += pow(Z[5], 2) * (6) * (-4) / den2_2;

            temp +=
                pow(Z[6], 2) * (8) * (108 * y2 - 36 * x2) * (1 + e2) / den2_3;
            temp +=
                pow(Z[7], 2) * (8) * (108 * x2 - 36 * y2) * (1 + e2) / den2_3;

            temp += pow(Z[8], 2) * (8) * (-36 * r2) / den2_4;
            temp += pow(Z[9], 2) * (8) * (-36 * r2) / den2_4;

            temp += pow(Z[10], 2) * (5) * 144 * (1 + e2 - 2 * r2) *
                    (1 + e2 - 6 * r2) / den2_5;

            temp += pow(Z[11], 2) * (10) * 36 *
                    (8 * (1 + e2 + e4) * x2 - 1 - e2 - e4 - e6) *
                    (1 + e2 + e4 + e6 - 8 * (1 + e2 + e4) * y2) * num2_6 /
                    den2_6;
            temp +=
                pow(Z[12], 2) * (10) * 36 *
                (-4 * pow(x_c - y_c, 2) * (e4 + e2 + 1) + 1 + e2 + e4 + e6) *
                (4 * pow(x_c + y_c, 2) * (e4 + e2 + 1) - 1 - e2 - e4 - e6) *
                num2_6 / den2_6;

            temp += pow(Z[13], 2) * (10) * (-144) * pow(r2, 2) / den2_7;
            temp += pow(Z[14], 2) * (10) * (-144) * pow(r2, 2) / den2_7;

            temp += pow(Z[15], 2) * (12) * 64 *
                    ((3 * e6 - 5 * e4 * r2 + 12 * e4 - 20 * e2 * r2 + 12 * e2 -
                      5 * r2 + 3) *
                     (9 * e6 * x2 - 3 * e6 * y2 - 25 * e4 * x4 -
                      20 * e4 * x2 * y2 + 36 * e4 * x2 + 5 * e4 * y4 -
                      12 * e4 * y2 - 100 * e2 * x4 - 80 * e2 * x2 * y2 +
                      36 * e2 * x2 + 20 * e2 * y4 - 12 * e2 * y2 - 25 * x4 -
                      20 * x2 * y2 + 9 * x2 + 5 * y4 - 3 * y2)) *
                    num2_8 / den2_8;
            temp += pow(Z[16], 2) * (12) * 64 *
                    (-(3 * e6 - 5 * e4 * r2 + 12 * e4 - 20 * e2 * r2 + 12 * e2 -
                       5 * r2 + 3) *
                     (3 * e6 * x2 - 9 * e6 * y2 - 5 * e4 * x4 +
                      20 * e4 * x2 * y2 + 12 * e4 * x2 + 25 * e4 * y4 -
                      36 * e4 * y2 - 20 * e2 * x4 + 80 * e2 * x2 * y2 +
                      12 * e2 * x2 + 100 * e2 * y4 - 36 * e2 * y2 - 5 * x4 +
                      20 * x2 * y2 + 3 * x2 + 25 * y4 - 9 * y2)) *
                    num2_8 / den2_8;

            temp +=
                pow(Z[17], 2) * (12) * 16.0 *
                (-36 * e16 * x2 - 36 * e16 * y2 + 180 * e14 * x4 +
                 360 * e14 * x2 * y2 - 72 * e14 * x2 + 180 * e14 * y4 -
                 72 * e14 * y2 - 125 * e12 * x6 - 1275 * e12 * x4 * y2 +
                 360 * e12 * x4 + 225 * e12 * x2 * y4 + 720 * e12 * x2 * y2 -
                 108 * e12 * x2 - 225 * e12 * y6 + 360 * e12 * y4 -
                 108 * e12 * y2 - 250 * e10 * x6 - 2550 * e10 * x4 * y2 +
                 540 * e10 * x4 + 450 * e10 * x2 * y4 + 1080 * e10 * x2 * y2 -
                 144 * e10 * x2 - 450 * e10 * y6 + 540 * e10 * y4 -
                 144 * e10 * y2 - 375 * e8 * x6 - 3825 * e8 * x4 * y2 +
                 720 * e8 * x4 + 675 * e8 * x2 * y4 + 1440 * e8 * x2 * y2 -
                 180 * e8 * x2 - 675 * e8 * y6 + 720 * e8 * y4 - 180 * e8 * y2 -
                 500 * e6 * x6 - 5100 * e6 * x4 * y2 + 720 * e6 * x4 +
                 900 * e6 * x2 * y4 + 1440 * e6 * x2 * y2 - 144 * e6 * x2 -
                 900 * e6 * y6 + 720 * e6 * y4 - 144 * e6 * y2 - 375 * e4 * x6 -
                 3825 * e4 * x4 * y2 + 540 * e4 * x4 + 675 * e4 * x2 * y4 +
                 1080 * e4 * x2 * y2 - 108 * e4 * x2 - 675 * e4 * y6 +
                 540 * e4 * y4 - 108 * e4 * y2 - 250 * e2 * x6 -
                 2550 * e2 * x4 * y2 + 360 * e2 * x4 + 450 * e2 * x2 * y4 +
                 720 * e2 * x2 * y2 - 72 * e2 * x2 - 450 * e2 * y6 +
                 360 * e2 * y4 - 72 * e2 * y2 - 125 * x6 - 1275 * x4 * y2 +
                 180 * x4 + 225 * x2 * y4 + 360 * x2 * y2 - 36 * x2 - 225 * y6 +
                 180 * y4 - 36 * y2) *
                num2_9 / den2_9;
            temp +=
                pow(Z[18], 2) * (12) * 16.0 *
                ((-225 * e12 - 450 * e10 - 675 * e8 - 900 * e6 - 675 * e4 -
                  450 * e2 - 225) *
                     x6 +
                 (180 * e14 + 225 * e12 * y2 + 360 * e12 + 450 * e10 * y2 +
                  540 * e10 + 675 * e8 * y2 + 720 * e8 + 900 * e6 * y2 +
                  720 * e6 + 675 * e4 * y2 + 540 * e4 + 450 * e2 * y2 +
                  360 * e2 + 225 * y2 + 180) *
                     x4 +
                 (-36 * e16 + 360 * e14 * y2 - 72 * e14 - 1275 * e12 * y4 +
                  720 * e12 * y2 - 108 * e12 - 2550 * e10 * y4 +
                  1080 * e10 * y2 - 144 * e10 - 3825 * e8 * y4 +
                  1440 * e8 * y2 - 180 * e8 - 5100 * e6 * y4 + 1440 * e6 * y2 -
                  144 * e6 - 3825 * e4 * y4 + 1080 * e4 * y2 - 108 * e4 -
                  2550 * e2 * y4 + 720 * e2 * y2 - 72 * e2 - 1275 * y4 +
                  360 * y2 - 36) *
                     x2 -
                 36 * e16 * y2 + 180 * e14 * y4 - 72 * e14 * y2 -
                 125 * e12 * y6 + 360 * e12 * y4 - 108 * e12 * y2 -
                 250 * e10 * y6 + 540 * e10 * y4 - 144 * e10 * y2 -
                 375 * e8 * y6 + 720 * e8 * y4 - 180 * e8 * y2 - 500 * e6 * y6 +
                 720 * e6 * y4 - 144 * e6 * y2 - 375 * e4 * y6 + 540 * e4 * y4 -
                 108 * e4 * y2 - 250 * e2 * y6 + 360 * e2 * y4 - 72 * e2 * y2 -
                 125 * y6 + 180 * y4 - 36 * y2) *
                num2_9 / den2_9;

            temp += pow(Z[19], 2) * (12) * (-400) * pow(r2, 3) / den2_10;
            temp += pow(Z[20], 2) * (12) * (-400) * pow(r2, 3) / den2_10;

            temp += pow(Z[21], 2) * (7) * 576 *
                    ((e4 - 5 * e2 * x2 - 5 * e2 * y2 + 3 * e2 + 5 * x4 +
                      10 * x2 * y2 - 5 * x2 + 5 * y4 - 5 * y2 + 1) *
                     (e4 - 15 * e2 * x2 - 15 * e2 * y2 + 3 * e2 + 25 * x4 +
                      50 * x2 * y2 - 15 * x2 + 25 * y4 - 15 * y2 + 1)) /
                    den2_11;

            out[ii] = temp;
        }
    }
    return result;
}

py::array_t<double> zernikeAnnularGrad(py::array_t<double> arrayZk,
                                       py::array_t<double> arrayX,
                                       py::array_t<double> arrayY, double e,
                                       std::string axis) {
    // Get the numpy array information
    ArrayNpInfo infoZk = getNpArrayInfo(arrayZk);
    ArrayNpInfo infoX = getNpArrayInfo(arrayX);
    ArrayNpInfo infoY = getNpArrayInfo(arrayY);

    // No pointer is passed, so NumPy will allocate the buffer
    size_t n = infoX.size;
    auto result = py::array_t<double>(n);
    ArrayNpInfo infoResult = getNpArrayInfo(result);

    // Assign the variables
    double *Z = infoZk.buf;
    double *x = infoX.buf;
    double *y = infoY.buf;
    double *d = infoResult.buf;

    // Parameters of constant
    double e2 = pow(e, 2);
    double e4 = e2 * e2;
    double e6 = e4 * e2;
    double e8 = e6 * e2;
    double e10 = e8 * e2;
    double e12 = e10 * e2;

    double sqrt_3 = sqrt(3);
    double sqrt_5 = sqrt(5);
    double sqrt_6 = sqrt(6);
    double sqrt_7 = sqrt(7);
    double sqrt_8 = sqrt(8);
    double sqrt_10 = sqrt(10);
    double sqrt_12 = sqrt(12);

    double den1 = sqrt(1 + e2);
    double den2 = 1 - e2;
    double den3 = sqrt(1 + e2 + e4);
    double den4 = sqrt(pow(1 - e2, 2) * (1 + e2) * (1 + 4 * e2 + e4));
    double den5 = sqrt(1 + e2 + e4 + e6);
    double den6 = pow(1 - e2, 2);

    double den7 = pow(1 - e2, 3) * (1 + e2 + e4);
    double num7 = sqrt(pow(1 - e2, 4) * (1 + e2 + e4) /
                       (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8));

    double den8 = sqrt(1 + e2 + e4 + e6 + e8);

    double den9 = pow(1 - e2, 3) * (1 + 4 * e2 + e4);
    double num9 =
        sqrt(pow(1 - e2, 2) * (1 + 4 * e2 + e4) / (1 + 9 * e2 + 9 * e4 + e6));

    double den10 = pow(1 - e2, 4) * (1 + e2) * (1 + e4);
    double num10 =
        sqrt(pow(1 - e2, 6) * (1 + e2) * (1 + e4) /
             (1 + 4 * e2 + 10 * e4 + 20 * e6 + 10 * e8 + 4 * e10 + e12));

    double den11 = sqrt(1 + e2 + e4 + e6 + e8 + e10);
    double den12 = pow(1 - e2, 3);

    // Parameters in loop
    double x2, y2, x4, y4, xy, r2, r4, temp, x_c, y_c;
    if (axis == "dx") {
        for (size_t ii = 0; ii < n; ii++) {
            x_c = x[ii];
            y_c = y[ii];

            x2 = x_c * x_c;
            y2 = y_c * y_c;
            x4 = x2 * x2;
            y4 = y2 * y2;
            xy = x_c * y_c;
            r2 = x2 + y2;

            // to make d an array with the same size as x
            temp = Z[0] * 0 * x_c;

            temp += Z[1] * 2 * 1 / den1;
            temp += Z[2] * 2 * 0;

            temp += Z[3] * sqrt_3 * 4 * x_c / den2;

            temp += Z[4] * sqrt_6 * 2 * y_c / den3;
            temp += Z[5] * sqrt_6 * 2 * x_c / den3;

            temp += Z[6] * sqrt_8 * 6 * xy * (1 + e2) / den4;
            temp += Z[7] * sqrt_8 *
                    ((9 * x2 + 3 * y2 - 2) * (1 + e2) - 2 * e4) / den4;

            temp += Z[8] * sqrt_8 * 6 * xy / den5;
            temp += Z[9] * sqrt_8 * (3 * x2 - 3 * y2) / den5;

            temp += Z[10] * sqrt_5 * 12 * x_c * (2 * r2 - 1 - e2) / den6;

            temp += Z[11] * sqrt_10 *
                    (x_c * (16 * x2 - 6) * (1 + e2 + e4) - 6 * x_c * e6) *
                    num7 / den7;
            temp +=
                Z[12] * sqrt_10 *
                (y_c * (24 * x2 + 8 * y2 - 6) * (1 + e2 + e4) - 6 * y_c * e6) *
                num7 / den7;

            temp += Z[13] * sqrt_10 * 4 * x_c * (x2 - 3 * y2) / den8;
            temp += Z[14] * sqrt_10 * 4 * y_c * (3 * x2 - y2) / den8;

            temp +=
                Z[15] * sqrt_12 *
                (3 * e8 - 36 * e6 * x2 - 12 * e6 * y2 + 12 * e6 + 50 * e4 * x4 +
                 60 * e4 * x2 * y2 - 144 * e4 * x2 + 10 * e4 * y4 -
                 48 * e4 * y2 + 30 * e4 + 200 * e2 * x4 + 240 * e2 * x2 * y2 -
                 144 * e2 * x2 + 40 * e2 * y4 - 48 * e2 * y2 + 12 * e2 +
                 50 * x4 + 60 * x2 * y2 - 36 * x2 + 10 * y4 - 12 * y2 + 3) *
                num9 / den9;
            temp += Z[16] * sqrt_12 *
                    (8 * xy *
                     (5 * r2 * (1 + 4 * e2 + e4) -
                      (3 + 12 * e2 + 12 * e4 + 3 * e6))) *
                    num9 / den9;

            temp += Z[17] * sqrt_12 *
                    (25 * (e6 + e4 + e2 + 1) * x4 +
                     (-12 * e8 - 30 * e6 * y2 - 12 * e6 - 30 * e4 * y2 -
                      12 * e4 - 30 * e2 * y2 - 12 * e2 - 30 * y2 - 12) *
                         x2 +
                     12 * e8 * y2 - 15 * e6 * y4 + 12 * e6 * y2 - 15 * e4 * y4 +
                     12 * e4 * y2 - 15 * e2 * y4 + 12 * e2 * y2 - 15 * y4 +
                     12 * y2) *
                    num10 / den10;
            temp +=
                Z[18] * sqrt_12 *
                (4.0 * xy *
                 (15 * (e6 + e4 + e2 + 1) * x2 - 6 * e8 + 5 * e6 * y2 - 6 * e6 +
                  5 * e4 * y2 - 6 * e4 + 5 * e2 * y2 - 6 * e2 + 5 * y2 - 6)) *
                num10 / den10;

            temp += Z[19] * sqrt_12 * 5 * (x2 * (x2 - 6 * y2) + y4) / den11;
            temp += Z[20] * sqrt_12 * 20 * xy * (x2 - y2) / den11;

            temp += Z[21] * sqrt_7 * 24 * x_c *
                    (e4 - e2 * (5 * y2 - 3) + 5 * x4 - 5 * y2 + 5 * y4 -
                     x2 * (5 * e2 - 10 * y2 + 5) + 1) /
                    den12;

            d[ii] = temp;
        }
    } else if (axis == "dy") {
        for (size_t ii = 0; ii < n; ii++) {
            x_c = x[ii];
            y_c = y[ii];

            x2 = x_c * x_c;
            y2 = y_c * y_c;
            x4 = x2 * x2;
            y4 = y2 * y2;
            xy = x_c * y_c;
            r2 = x2 + y2;

            // to make d an array with the same size as x
            temp = Z[0] * 0 * x_c;

            temp += Z[1] * 2 * 0;
            temp += Z[2] * 2 * 1 / den1;

            temp += Z[3] * sqrt_3 * 4 * y_c / den2;

            temp += Z[4] * sqrt_6 * 2 * x_c / den3;
            temp += Z[5] * sqrt_6 * (-2) * y_c / den3;

            temp += Z[6] * sqrt_8 *
                    ((1 + e2) * (3 * x2 + 9 * y2 - 2) - 2 * e4) / den4;
            temp += Z[7] * sqrt_8 * 6 * xy * (1 + e2) / den4;

            temp += Z[8] * sqrt_8 * (3 * x2 - 3 * y2) / den5;
            temp += Z[9] * sqrt_8 * (-6) * xy / den5;

            temp += Z[10] * sqrt_5 * 12 * y_c * (2 * r2 - 1 - e2) / den6;

            temp += Z[11] * sqrt_10 *
                    (y_c * (6 - 16 * y2) * (1 + e2 + e4) + 6 * y_c * e6) *
                    num7 / den7;
            temp +=
                Z[12] * sqrt_10 *
                (x_c * (8 * x2 + 24 * y2 - 6) * (1 + e2 + e4) - 6 * x_c * e6) *
                num7 / den7;

            temp += Z[13] * sqrt_10 * 4 * y_c * (y2 - 3 * x2) / den8;
            temp += Z[14] * sqrt_10 * 4 * x_c * (x2 - 3 * y2) / den8;

            temp += Z[15] * sqrt_12 *
                    (-x_c * (24 * y_c + 4 * e2 * (24 * y_c - 40 * y_c * r2) +
                             2 * e4 * (48 * y_c - 20 * y_c * r2) +
                             24 * e6 * y_c - 40 * y_c * r2)) *
                    num9 / den9;
            temp +=
                Z[16] * sqrt_12 *
                (3 * e8 - 12 * e6 * x2 - 36 * e6 * y2 + 12 * e6 + 10 * e4 * x4 +
                 60 * e4 * x2 * y2 - 48 * e4 * x2 + 50 * e4 * y4 -
                 144 * e4 * y2 + 30 * e4 + 40 * e2 * x4 + 240 * e2 * x2 * y2 -
                 48 * e2 * x2 + 200 * e2 * y4 - 144 * e2 * y2 + 12 * e2 +
                 10 * x4 + 60 * x2 * y2 - 12 * x2 + 50 * y4 - 36 * y2 + 3) *
                num9 / den9;

            temp += Z[17] * sqrt_12 *
                    (4.0 * xy *
                     ((-5) * (e6 + e4 + e2 + 1) * x2 + 6 * e8 - 15 * e6 * y2 +
                      6 * e6 - 15 * e4 * y2 + 6 * e4 - 15 * e2 * y2 + 6 * e2 -
                      15 * y2 + 6)) *
                    num10 / den10;
            temp += Z[18] * sqrt_12 *
                    (-12 * e8 * x2 + 12 * e8 * y2 + 15 * e6 * x4 +
                     30 * e6 * x2 * y2 - 12 * e6 * x2 - 25 * e6 * y4 +
                     12 * e6 * y2 + 15 * e4 * x4 + 30 * e4 * x2 * y2 -
                     12 * e4 * x2 - 25 * e4 * y4 + 12 * e4 * y2 + 15 * e2 * x4 +
                     30 * e2 * x2 * y2 - 12 * e2 * x2 - 25 * e2 * y4 +
                     12 * e2 * y2 + 15 * x4 + 30 * x2 * y2 - 12 * x2 - 25 * y4 +
                     12 * y2) *
                    num10 / den10;

            temp += Z[19] * sqrt_12 * 20 * xy * (y2 - x2) / den11;
            temp += Z[20] * sqrt_12 * 5 * (x2 * (x2 - 6 * y2) + y4) / den11;

            temp += Z[21] * sqrt_7 * 24 * y_c *
                    (e4 - e2 * (5 * x2 - 3) - 5 * x2 + 5 * x4 + 5 * y4 -
                     y2 * (5 * e2 - 10 * x2 + 5) + 1) /
                    den12;

            d[ii] = temp;
        }
    } else if (axis == "dx2") {
        for (size_t ii = 0; ii < n; ii++) {
            x_c = x[ii];
            y_c = y[ii];

            x2 = x_c * x_c;
            y2 = y_c * y_c;
            x4 = x2 * x2;
            y4 = y2 * y2;
            xy = x_c * y_c;
            r2 = x2 + y2;
            r4 = r2 * r2;

            // to make d an array with the same size as x
            temp = Z[0] * 0 * x_c;
            temp += Z[1] * 0;
            temp += Z[2] * 0;

            temp += Z[3] * sqrt_3 * 4 / den2;
            temp += Z[4] * 0;

            temp += Z[5] * sqrt_6 * 2 / den3;

            temp += Z[6] * sqrt_8 * 6 * y_c * (1 + e2) / den4;
            temp += Z[7] * sqrt_8 * 18 * x_c * (1 + e2) / den4;

            temp += Z[8] * sqrt_8 * 6 * y_c / den5;
            temp += Z[9] * sqrt_8 * 6 * x_c / den5;

            temp += Z[10] * sqrt_5 * 12 * (6 * x2 + 2 * y2 - e2 - 1) / den6;

            temp += Z[11] * sqrt_10 * ((48 * x2 - 6) * (1 + e2 + e4) - 6 * e6) *
                    num7 / den7;
            temp += Z[12] * sqrt_10 * 48 * xy * (1 + e2 + e4) * num7 / den7;

            temp += Z[13] * sqrt_10 * 12 * (x2 - y2) / den8;
            temp += Z[14] * sqrt_10 * 24 * xy / den8;

            temp += Z[15] * sqrt_12 *
                    (-8 * x_c *
                     (9 * e6 - 25 * e4 * x2 - 15 * e4 * y2 + 36 * e4 -
                      100 * e2 * x2 - 60 * e2 * y2 + 36 * e2 - 25 * x2 -
                      15 * y2 + 9)) *
                    num9 / den9;
            temp +=
                Z[16] * sqrt_12 *
                (-8 * y_c *
                 (3 * e6 - 15 * e4 * x2 - 5 * e4 * y2 + 12 * e4 - 60 * e2 * x2 -
                  20 * e2 * y2 + 12 * e2 - 15 * x2 - 5 * y2 + 3)) *
                num9 / den9;

            temp += Z[17] * sqrt_12 *
                    (-4 * x_c *
                     (6 * e8 - 25 * e6 * x2 + 15 * e6 * y2 + 6 * e6 -
                      25 * e4 * x2 + 15 * e4 * y2 + 6 * e4 - 25 * e2 * x2 +
                      15 * e2 * y2 + 6 * e2 - 25 * x2 + 15 * y2 + 6)) *
                    num10 / den10;
            temp += Z[18] * sqrt_12 *
                    (-4 * y_c *
                     (6 * e8 - 45 * e6 * x2 - 5 * e6 * y2 + 6 * e6 -
                      45 * e4 * x2 - 5 * e4 * y2 + 6 * e4 - 45 * e2 * x2 -
                      5 * e2 * y2 + 6 * e2 - 45 * x2 - 5 * y2 + 6)) *
                    num10 / den10;

            temp += Z[19] * sqrt_12 * 20 * x_c * (x2 - 3 * y2) / den11;
            temp += Z[20] * sqrt_12 * 20 * y_c * (3 * x2 - y2) / den11;

            temp += Z[21] * sqrt_7 *
                    (480 * x2 * r2 + 120 * r4 + 24 * e4 - 360 * x2 - 120 * y2 -
                     3 * e2 * (120 * x2 + 40 * y2 - 24) + 24) /
                    den12;

            d[ii] = temp;
        }
    } else if (axis == "dy2") {
        for (size_t ii = 0; ii < n; ii++) {
            x_c = x[ii];
            y_c = y[ii];

            x2 = x_c * x_c;
            y2 = y_c * y_c;
            x4 = x2 * x2;
            y4 = y2 * y2;
            xy = x_c * y_c;
            r2 = x2 + y2;
            r4 = r2 * r2;

            // to make d an array with the same size as x
            temp = Z[0] * 0 * x_c;
            temp += Z[1] * 0;
            temp += Z[2] * 0;

            temp += Z[3] * sqrt_3 * 4 / den2;
            temp += Z[4] * 0;

            temp += Z[5] * sqrt_6 * (-2) / den3;

            temp += Z[6] * sqrt_8 * (1 + e2) * 18 * y_c / den4;
            temp += Z[7] * sqrt_8 * 6 * x_c * (1 + e2) / den4;

            temp += Z[8] * sqrt_8 * (-6) * y_c / den5;
            temp += Z[9] * sqrt_8 * (-6) * x_c / den5;

            temp += Z[10] * sqrt_5 * 12 * (2 * x2 + 6 * y2 - e2 - 1) / den6;

            temp += Z[11] * sqrt_10 * ((6 - 48 * y2) * (1 + e2 + e4) + 6 * e6) *
                    num7 / den7;
            temp += Z[12] * sqrt_10 * 48 * xy * (1 + e2 + e4) * num7 / den7;

            temp += Z[13] * sqrt_10 * 12 * (y2 - x2) / den8;
            temp += Z[14] * sqrt_10 * (-24) * xy / den8;

            temp +=
                Z[15] * sqrt_12 *
                (-8 * x_c *
                 (3 * e6 - 5 * e4 * x2 - 15 * e4 * y2 + 12 * e4 - 20 * e2 * x2 -
                  60 * e2 * y2 + 12 * e2 - 5 * x2 - 15 * y2 + 3)) *
                num9 / den9;
            temp += Z[16] * sqrt_12 *
                    (-8 * y_c *
                     (9 * e6 - 15 * e4 * x2 - 25 * e4 * y2 + 36 * e4 -
                      60 * e2 * x2 - 100 * e2 * y2 + 36 * e2 - 15 * x2 -
                      25 * y2 + 9)) *
                    num9 / den9;

            temp += Z[17] * sqrt_12 *
                    (4 * x_c *
                     (6 * e8 - 5 * e6 * x2 - 45 * e6 * y2 + 6 * e6 -
                      5 * e4 * x2 - 45 * e4 * y2 + 6 * e4 - 5 * e2 * x2 -
                      45 * e2 * y2 + 6 * e2 - 5 * x2 - 45 * y2 + 6)) *
                    num10 / den10;
            temp += Z[18] * sqrt_12 *
                    (4 * y_c *
                     (6 * e8 + 15 * e6 * x2 - 25 * e6 * y2 + 6 * e6 +
                      15 * e4 * x2 - 25 * e4 * y2 + 6 * e4 + 15 * e2 * x2 -
                      25 * e2 * y2 + 6 * e2 + 15 * x2 - 25 * y2 + 6)) *
                    num10 / den10;

            temp += Z[19] * sqrt_12 * 20 * x_c * (3 * y2 - x2) / den11;
            temp += Z[20] * sqrt_12 * 20 * y_c * (y2 - 3 * x2) / den11;

            temp += Z[21] * sqrt_7 *
                    (480 * y2 * r2 + 120 * r4 + 24 * e4 - 120 * x2 - 360 * y2 -
                     3 * e2 * (40 * x2 + 120 * y2 - 24) + 24) /
                    den12;

            d[ii] = temp;
        }
    } else if (axis == "dxy") {
        for (size_t ii = 0; ii < n; ii++) {
            x_c = x[ii];
            y_c = y[ii];

            x2 = x_c * x_c;
            y2 = y_c * y_c;
            x4 = x2 * x2;
            y4 = y2 * y2;
            xy = x_c * y_c;
            r2 = x2 + y2;
            r4 = r2 * r2;

            // to make d an array with the same size as x
            temp = Z[0] * 0 * x_c;
            temp += Z[1] * 0;
            temp += Z[2] * 0;
            temp += Z[3] * 0;

            temp += Z[4] * sqrt_6 * 2 / den3;
            temp += Z[5] * 0;

            temp += Z[6] * sqrt_8 * (1 + e2) * (6 * x_c) / den4;
            temp += Z[7] * sqrt_8 * 6 * y_c * (1 + e2) / den4;

            temp += Z[8] * sqrt_8 * 6 * x_c / den5;
            temp += Z[9] * sqrt_8 * (-6) * y_c / den5;

            temp += Z[10] * sqrt_5 * 48 * xy / den6;

            temp += Z[11] * sqrt_10 * 0;
            temp += Z[12] * sqrt_10 *
                    ((24 * x2 + 24 * y2 - 6) * (1 + e2 + e4) - 6 * e6) * num7 /
                    den7;

            temp += Z[13] * sqrt_10 * (-24) * xy / den8;
            temp += Z[14] * sqrt_10 * 12 * (x2 - y2) / den8;

            temp +=
                Z[15] * sqrt_12 *
                (-8 * y_c *
                 (3 * e6 - 15 * e4 * x2 - 5 * e4 * y2 + 12 * e4 - 60 * e2 * x2 -
                  20 * e2 * y2 + 12 * e2 - 15 * x2 - 5 * y2 + 3)) *
                num9 / den9;
            temp +=
                Z[16] * sqrt_12 *
                (-8 * x_c *
                 (3 * e6 - 5 * e4 * x2 - 15 * e4 * y2 + 12 * e4 - 20 * e2 * x2 -
                  60 * e2 * y2 + 12 * e2 - 5 * x2 - 15 * y2 + 3)) *
                num9 / den9;

            temp += Z[17] * sqrt_12 *
                    (12 * y_c *
                     (2 * e8 - 5 * e6 * r2 + 2 * e6 - 5 * e4 * r2 + 2 * e4 -
                      5 * e2 * r2 + 2 * e2 - 5 * r2 + 2)) *
                    num10 / den10;
            temp += Z[18] * sqrt_12 *
                    (-12 * x_c *
                     (2 * e8 - 5 * e6 * r2 + 2 * e6 - 5 * e4 * r2 + 2 * e4 -
                      5 * e2 * r2 + 2 * e2 - 5 * r2 + 2)) *
                    num10 / den10;

            temp += Z[19] * sqrt_12 * 20 * y_c * (y2 - 3 * x2) / den11;
            temp += Z[20] * sqrt_12 * 20 * x_c * (x2 - 3 * y2) / den11;

            temp += Z[21] * sqrt_7 * 240 * xy * (2 * r2 - 1 - e2) / den12;

            d[ii] = temp;
        }
    }
    return result;
}

py::array_t<double> poly10_2D(py::array_t<double> arrayC,
                              py::array_t<double> arrayX,
                              py::array_t<double> arrayY) {
    // Get the numpy array information
    ArrayNpInfo infoC = getNpArrayInfo(arrayC);
    ArrayNpInfo infoX = getNpArrayInfo(arrayX);
    ArrayNpInfo infoY = getNpArrayInfo(arrayY);

    // No pointer is passed, so NumPy will allocate the buffer
    size_t n = infoX.size;
    auto result = py::array_t<double>(n);
    ArrayNpInfo infoResult = getNpArrayInfo(result);

    // Assign the variables
    double *c = infoC.buf;
    double *x = infoX.buf;
    double *y = infoY.buf;
    double *cyOut = infoResult.buf;

    double x_c, y_c;
    for (size_t ii = 0; ii < n; ii++) {
        x_c = x[ii];
        y_c = y[ii];
        cyOut[ii] =
            c[0] + c[1] * x_c + c[2] * y_c + c[3] * x_c * x_c +
            c[4] * x_c * y_c + c[5] * y_c * y_c + c[6] * pow(x_c, 3) +
            c[7] * pow(x_c, 2) * y_c + c[8] * x_c * pow(y_c, 2) +
            c[9] * pow(y_c, 3) + c[10] * pow(x_c, 4) +
            c[11] * pow(x_c, 3) * y_c + c[12] * pow(x_c, 2) * pow(y_c, 2) +
            c[13] * x_c * pow(y_c, 3) + c[14] * pow(y_c, 4) +
            c[15] * pow(x_c, 5) + c[16] * pow(x_c, 4) * y_c +
            c[17] * pow(x_c, 3) * pow(y_c, 2) +
            c[18] * pow(x_c, 2) * pow(y_c, 3) + c[19] * x_c * pow(y_c, 4) +
            c[20] * pow(y_c, 5) + c[21] * pow(x_c, 6) +
            c[22] * pow(x_c, 5) * y_c + c[23] * pow(x_c, 4) * pow(y_c, 2) +
            c[24] * pow(x_c, 3) * pow(y_c, 3) +
            c[25] * pow(x_c, 2) * pow(y_c, 4) + c[26] * x_c * pow(y_c, 5) +
            c[27] * pow(y_c, 6) + c[28] * pow(x_c, 7) +
            c[29] * pow(x_c, 6) * y_c + c[30] * pow(x_c, 5) * pow(y_c, 2) +
            c[31] * pow(x_c, 4) * pow(y_c, 3) +
            c[32] * pow(x_c, 3) * pow(y_c, 4) +
            c[33] * pow(x_c, 2) * pow(y_c, 5) + c[34] * x_c * pow(y_c, 6) +
            c[35] * pow(y_c, 7) + c[36] * pow(x_c, 8) +
            c[37] * pow(x_c, 7) * y_c + c[38] * pow(x_c, 6) * pow(y_c, 2) +
            c[39] * pow(x_c, 5) * pow(y_c, 3) +
            c[40] * pow(x_c, 4) * pow(y_c, 4) +
            c[41] * pow(x_c, 3) * pow(y_c, 5) +
            c[42] * pow(x_c, 2) * pow(y_c, 6) + c[43] * x_c * pow(y_c, 7) +
            c[44] * pow(y_c, 8) + c[45] * pow(x_c, 9) +
            c[46] * pow(x_c, 8) * y_c + c[47] * pow(x_c, 7) * pow(y_c, 2) +
            c[48] * pow(x_c, 6) * pow(y_c, 3) +
            c[49] * pow(x_c, 5) * pow(y_c, 4) +
            c[50] * pow(x_c, 4) * pow(y_c, 5) +
            c[51] * pow(x_c, 3) * pow(y_c, 6) +
            c[52] * pow(x_c, 2) * pow(y_c, 7) + c[53] * x_c * pow(y_c, 8) +
            c[54] * pow(y_c, 9) + c[55] * pow(x_c, 10) +
            c[56] * pow(x_c, 9) * y_c + c[57] * pow(x_c, 8) * pow(y_c, 2) +
            c[58] * pow(x_c, 7) * pow(y_c, 3) +
            c[59] * pow(x_c, 6) * pow(y_c, 4) +
            c[60] * pow(x_c, 5) * pow(y_c, 5) +
            c[61] * pow(x_c, 4) * pow(y_c, 6) +
            c[62] * pow(x_c, 3) * pow(y_c, 7) +
            c[63] * pow(x_c, 2) * pow(y_c, 8) + c[64] * x_c * pow(y_c, 9) +
            c[65] * pow(y_c, 10);
    }
    return result;
}

py::array_t<double> poly10Grad(py::array_t<double> arrayC,
                               py::array_t<double> arrayX,
                               py::array_t<double> arrayY, std::string axis) {
    // Get the numpy array information
    ArrayNpInfo infoC = getNpArrayInfo(arrayC);
    ArrayNpInfo infoX = getNpArrayInfo(arrayX);
    ArrayNpInfo infoY = getNpArrayInfo(arrayY);

    // No pointer is passed, so NumPy will allocate the buffer
    size_t n = infoX.size;
    auto result = py::array_t<double>(n);
    ArrayNpInfo infoResult = getNpArrayInfo(result);

    // Assign the variables
    double *c = infoC.buf;
    double *x = infoX.buf;
    double *y = infoY.buf;
    double *cy_out = infoResult.buf;

    double x_c, y_c;
    if (axis == "dx") {
        for (size_t ii = 0; ii < n; ii++) {
            x_c = x[ii];
            y_c = y[ii];
            cy_out[ii] =
                c[1] + c[3] * 2 * x_c + c[4] * y_c + c[6] * 3 * pow(x_c, 2) +
                c[7] * 2 * x_c * y_c + c[8] * pow(y_c, 2) +
                c[10] * 4 * pow(x_c, 3) + c[11] * 3 * pow(x_c, 2) * y_c +
                c[12] * 2 * x_c * pow(y_c, 2) + c[13] * pow(y_c, 3) +
                c[15] * 5 * pow(x_c, 4) + c[16] * 4 * pow(x_c, 3) * y_c +
                c[17] * 3 * pow(x_c, 2) * pow(y_c, 2) +
                c[18] * 2 * x_c * pow(y_c, 3) + c[19] * pow(y_c, 4) +
                c[21] * 6 * pow(x_c, 5) + c[22] * 5 * pow(x_c, 4) * y_c +
                c[23] * 4 * pow(x_c, 3) * pow(y_c, 2) +
                c[24] * 3 * pow(x_c, 2) * pow(y_c, 3) +
                c[25] * 2 * x_c * pow(y_c, 4) + c[26] * pow(y_c, 5) +
                c[28] * 7 * pow(x_c, 6) + c[29] * 6 * pow(x_c, 5) * y_c +
                c[30] * 5 * pow(x_c, 4) * pow(y_c, 2) +
                c[31] * 4 * pow(x_c, 3) * pow(y_c, 3) +
                c[32] * 3 * pow(x_c, 2) * pow(y_c, 4) +
                c[33] * 2 * x_c * pow(y_c, 5) + c[34] * pow(y_c, 6) +
                c[36] * 8 * pow(x_c, 7) + c[37] * 7 * pow(x_c, 6) * y_c +
                c[38] * 6 * pow(x_c, 5) * pow(y_c, 2) +
                c[39] * 5 * pow(x_c, 4) * pow(y_c, 3) +
                c[40] * 4 * pow(x_c, 3) * pow(y_c, 4) +
                c[41] * 3 * pow(x_c, 2) * pow(y_c, 5) +
                c[42] * 2 * x_c * pow(y_c, 6) + c[43] * pow(y_c, 7) +
                c[45] * 9 * pow(x_c, 8) + c[46] * 8 * pow(x_c, 7) * y_c +
                c[47] * 7 * pow(x_c, 6) * pow(y_c, 2) +
                c[48] * 6 * pow(x_c, 5) * pow(y_c, 3) +
                c[49] * 5 * pow(x_c, 4) * pow(y_c, 4) +
                c[50] * 4 * pow(x_c, 3) * pow(y_c, 5) +
                c[51] * 3 * pow(x_c, 2) * pow(y_c, 6) +
                c[52] * 2 * x_c * pow(y_c, 7) + c[53] * pow(y_c, 8) +
                c[55] * 10 * pow(x_c, 9) + c[56] * 9 * pow(x_c, 8) * y_c +
                c[57] * 8 * pow(x_c, 7) * pow(y_c, 2) +
                c[58] * 7 * pow(x_c, 6) * pow(y_c, 3) +
                c[59] * 6 * pow(x_c, 5) * pow(y_c, 4) +
                c[60] * 5 * pow(x_c, 4) * pow(y_c, 5) +
                c[61] * 4 * pow(x_c, 3) * pow(y_c, 6) +
                c[62] * 3 * pow(x_c, 2) * pow(y_c, 7) +
                c[63] * 2 * x_c * pow(y_c, 8) + c[64] * pow(y_c, 9);
        }
    } else if (axis == "dy") {
        for (size_t ii = 0; ii < n; ii++) {
            x_c = x[ii];
            y_c = y[ii];
            cy_out[ii] =
                c[2] + c[4] * x_c + c[5] * 2 * y_c + c[7] * pow(x_c, 2) +
                c[8] * x_c * 2 * y_c + c[9] * 3 * pow(y_c, 2) +
                c[11] * pow(x_c, 3) + c[12] * pow(x_c, 2) * 2 * y_c +
                c[13] * x_c * 3 * pow(y_c, 2) + c[14] * 4 * pow(y_c, 3) +
                c[16] * pow(x_c, 4) + c[17] * pow(x_c, 3) * 2 * y_c +
                c[18] * pow(x_c, 2) * 3 * pow(y_c, 2) +
                c[19] * x_c * 4 * pow(y_c, 3) + c[20] * 5 * pow(y_c, 4) +
                c[22] * pow(x_c, 5) + c[23] * pow(x_c, 4) * 2 * y_c +
                c[24] * pow(x_c, 3) * 3 * pow(y_c, 2) +
                c[25] * pow(x_c, 2) * 4 * pow(y_c, 3) +
                c[26] * x_c * 5 * pow(y_c, 4) + c[27] * 6 * pow(y_c, 5) +
                c[29] * pow(x_c, 6) + c[30] * pow(x_c, 5) * 2 * y_c +
                c[31] * pow(x_c, 4) * 3 * pow(y_c, 2) +
                c[32] * pow(x_c, 3) * 4 * pow(y_c, 3) +
                c[33] * pow(x_c, 2) * 5 * pow(y_c, 4) +
                c[34] * x_c * 6 * pow(y_c, 5) + c[35] * 7 * pow(y_c, 6) +
                c[37] * pow(x_c, 7) + c[38] * pow(x_c, 6) * 2 * y_c +
                c[39] * pow(x_c, 5) * 3 * pow(y_c, 2) +
                c[40] * pow(x_c, 4) * 4 * pow(y_c, 3) +
                c[41] * pow(x_c, 3) * 5 * pow(y_c, 4) +
                c[42] * pow(x_c, 2) * 6 * pow(y_c, 5) +
                c[43] * x_c * 7 * pow(y_c, 6) + c[44] * 8 * pow(y_c, 7) +
                c[46] * pow(x_c, 8) + c[47] * pow(x_c, 7) * 2 * y_c +
                c[48] * pow(x_c, 6) * 3 * pow(y_c, 2) +
                c[49] * pow(x_c, 5) * 4 * pow(y_c, 3) +
                c[50] * pow(x_c, 4) * 5 * pow(y_c, 4) +
                c[51] * pow(x_c, 3) * 6 * pow(y_c, 5) +
                c[52] * pow(x_c, 2) * 7 * pow(y_c, 6) +
                c[53] * x_c * 8 * pow(y_c, 7) + c[54] * 9 * pow(y_c, 8) +
                c[56] * pow(x_c, 9) + c[57] * pow(x_c, 8) * 2 * y_c +
                c[58] * pow(x_c, 7) * 3 * pow(y_c, 2) +
                c[59] * pow(x_c, 6) * 4 * pow(y_c, 3) +
                c[60] * pow(x_c, 5) * 5 * pow(y_c, 4) +
                c[61] * pow(x_c, 4) * 6 * pow(y_c, 5) +
                c[62] * pow(x_c, 3) * 7 * pow(y_c, 6) +
                c[63] * pow(x_c, 2) * 8 * pow(y_c, 7) +
                c[64] * x_c * 9 * pow(y_c, 8) + c[65] * 10 * pow(y_c, 9);
        }
    }
    return result;
}

} // namespace cwfs
} // namespace wep
} // namespace ts
} // namespace lsst
