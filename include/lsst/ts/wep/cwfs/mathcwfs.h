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

#ifndef _LSST_TS_WEP_CWFS_MATHCWFS_H
#define _LSST_TS_WEP_CWFS_MATHCWFS_H

#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace lsst {
namespace ts {
namespace wep {
namespace cwfs {

/**
 * Annular Zernike polynomials evaluation.
 *
 * @param[in] arrayZk  Coefficient of annular Zernike polynomials
 * @param[in] arrayX  X coordinate on pupil plane
 * @param[in] arrayY  Y coordinate on pupil plane
 * @param[in] e  Obscuration value
 * @return the wavefront surface
 */
py::array_t<double> zernikeAnnularEval(py::array_t<double> arrayZk,
                                       py::array_t<double> arrayX,
                                       py::array_t<double> arrayY, double e);

/**
 * Jacobian of annular Zernike polynomials.
 *
 * @param[in] arrayZk  Coefficient of annular Zernike polynomials
 * @param[in] arrayX  X coordinate on pupil plane
 * @param[in] arrayY  Y coordinate on pupil plane
 * @param[in] e  Obscuration value
 * @param[in] atype  Type/ Order of Jocobian Matrix ("1st" or "2nd")
 * @return the Jacobian elements in pupul x and y directions
 */
py::array_t<double> zernikeAnnularJacobian(py::array_t<double> arrayZk,
                                           py::array_t<double> arrayX,
                                           py::array_t<double> arrayY, double e,
                                           std::string atype);

/**
 * Gradient of annular Zernike polynomials.
 *
 * @param[in] arrayZk  Coefficient of annular Zernike polynomials
 * @param[in] arrayX  X coordinate on pupil plane
 * @param[in] arrayY  Y coordinate on pupil plane
 * @param[in] e  Obscuration value
 * @param[in] axis  Axis of "dx", "dy", "dx2", "dy2", or "dxy"
 * @return the integration elements of gradient in pupul x and y directions
 */
py::array_t<double> zernikeAnnularGrad(py::array_t<double> arrayZk,
                                       py::array_t<double> arrayX,
                                       py::array_t<double> arrayY, double e,
                                       std::string axis);

/**
 * Polynomial fit to 10th order in 2D (x, y dimensions).
 *
 * @param[in] arrayC  Parameters of off-axis distrotion
 * @param[in] arrayX  X coordinate on pupil plane
 * @param[in] arrayY  Y coordinate on pupil plane
 * @return the corrected parameters for off-axis distortion
 */
py::array_t<double> poly10_2D(py::array_t<double> arrayC,
                              py::array_t<double> arrayX,
                              py::array_t<double> arrayY);

/**
 * Gradient of polynomial fit to 10th order in 2D (x, y dimensions).
 *
 * @param[in] arrayC  Parameters of off-axis distrotion
 * @param[in] arrayX  X coordinate on pupil plane
 * @param[in] arrayY  Y coordinate on pupil plane
 * @param[in] axis  Direction of gradient ("dx" or "dy")
 * @return the corrected parameters for off-axis distortion
 */
py::array_t<double> poly10Grad(py::array_t<double> arrayC,
                               py::array_t<double> arrayX,
                               py::array_t<double> arrayY, std::string axis);

PYBIND11_MODULE(mathcwfs, m) {
    m.def("zernikeAnnularEval", &zernikeAnnularEval,
          "Jacobian of annular Zernike polynomials.");
    m.def("zernikeAnnularJacobian", &zernikeAnnularJacobian,
          "Jacobian of annular Zernike polynomials.");
    m.def("zernikeAnnularGrad", &zernikeAnnularGrad,
          "Gradient of annular Zernike polynomials.");
    m.def("poly10_2D", &poly10_2D,
          "Polynomial fit to 10th order in 2D (x, y dimensions).");
    m.def("poly10Grad", &poly10Grad,
          "Gradient of polynomial fit to 10th order in 2D (x, y dimensions).");
}

} // namespace cwfs
} // namespace wep
} // namespace ts
} // namespace lsst

#endif // !defined(_LSST_TS_WEP_CWFS_MATHCWFS_H)
