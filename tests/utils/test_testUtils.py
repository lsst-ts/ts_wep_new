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

import os
import unittest

import threadpoolctl
from lsst.ts.wep.utils.testUtils import enforce_single_threading


class TestSingleThreading(unittest.TestCase):
    def setUp(self):
        # Store original environment variables to restore after test
        self.original_openblas = os.environ.get("OPENBLAS_NUM_THREADS")
        self.original_mkl = os.environ.get("MKL_NUM_THREADS")
        self.original_omp = os.environ.get("OMP_NUM_THREADS")

    def test_enforce_single_threading(self):
        # Call the function to enforce single-threading
        enforce_single_threading()

        # Check environment variables
        self.assertEqual(
            os.environ["OPENBLAS_NUM_THREADS"],
            "1",
            "OPENBLAS_NUM_THREADS should be set to 1",
        )
        self.assertEqual(
            os.environ["MKL_NUM_THREADS"], "1", "MKL_NUM_THREADS should be set to 1"
        )
        self.assertEqual(
            os.environ["OMP_NUM_THREADS"], "1", "OMP_NUM_THREADS should be set to 1"
        )

    def test_threadpool_limits(self):
        # Call the function to enforce single-threading
        enforce_single_threading()

        # Get current threadpool limits
        threadpool_info = threadpoolctl.threadpool_info()

        # Check that BLAS libraries are limited to 1 thread
        for library in threadpool_info:
            if library["user_api"] == "blas":
                self.assertEqual(
                    library["num_threads"],
                    1,
                    f"BLAS library {library['prefix']} should be limited to 1 thread",
                )

    def tearDown(self):
        # Restore original environment variables
        if self.original_openblas is not None:
            os.environ["OPENBLAS_NUM_THREADS"] = self.original_openblas
        else:
            os.environ.pop("OPENBLAS_NUM_THREADS", None)

        if self.original_mkl is not None:
            os.environ["MKL_NUM_THREADS"] = self.original_mkl
        else:
            os.environ.pop("MKL_NUM_THREADS", None)

        if self.original_omp is not None:
            os.environ["OMP_NUM_THREADS"] = self.original_omp
        else:
            os.environ.pop("OMP_NUM_THREADS", None)


if __name__ == "__main__":
    unittest.main()
