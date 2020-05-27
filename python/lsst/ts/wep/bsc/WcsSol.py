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

import numpy as np

from lsst.sims.utils import ObservationMetaData
from lsst.obs.lsstSim import LsstSimMapper
from lsst.sims.coordUtils.CameraUtils import (
    raDecFromPixelCoords,
    pixelCoordsFromRaDec,
    focalPlaneCoordsFromRaDec,
)


class WcsSol(object):
    def __init__(self, camera=None):
        """Initialize the world coordinate system (WCS) solution class.

        Parameters
        ----------
        camera : lsst.afw.cameraGeom.camera.camera.Camera, optional
            A collection of Detectors that also supports coordinate
            transformation. (the default is None.)
        """

        self._obs = ObservationMetaData()

        if camera is None:
            self._camera = LsstSimMapper().camera
        else:
            self._camera = camera

    def setCamera(self, camera):
        """Set the camera object.

        Parameters
        ----------
        camera : lsst.afw.cameraGeom.camera.camera.Camera
            A collection of Detectors that also supports coordinate
            transformation.
        """

        self._camera = camera

    def getCamera(self):
        """Get the camera object.

        Returns
        -------
        lsst.afw.cameraGeom.camera.camera.Camera
            A collection of Detectors that also supports coordinate
            transformation.
        """

        return self._camera

    def setObsMetaData(self, ra, dec, rotSkyPos, mjd=59580.0):
        """Set the observation meta data.

        Parameters
        ----------
        ra : float
            Pointing ra in degree.
        dec : float
            Pointing decl in degree.
        rotSkyPos : float
            The orientation of the telescope in degrees.
        mjd : float
            Camera MJD. (the default is 59580.0.)
        """

        self._obs = ObservationMetaData(
            pointingRA=ra, pointingDec=dec, rotSkyPos=rotSkyPos, mjd=mjd
        )

    def raDecFromPixelCoords(
        self, xPix, yPix, chipName, epoch=2000.0, includeDistortion=True
    ):
        """Convert pixel coordinates into RA, Dec.

        WARNING: This method does not account for apparent motion due to
        parallax. This method is only useful for mapping positions on a
        theoretical focal plane to positions on the celestial sphere.

        Parameters
        ----------
        xPix : float or numpy.ndarray
            xPix is the x pixel coordinate.
        yPix : float or numpy.ndarray
            yPix is the y pixel coordinate.
        chipName : str or numpy.ndarray
            chipName is the name of the chip(s) on which the pixel coordinates
            are defined.  This can be an array (in which case there should be
            one chip name for each (xPix, yPix) coordinate pair), or a single
            value (in which case, all of the (xPix, yPix) points will be
            reckoned on that chip).
        epoch : float, optional
            epoch is the mean epoch in years of the celestial coordinate
            system. (the default is 2000.0.)
        includeDistortion : bool, optional
            If True (default), then this method will expect the true pixel
            coordinates with optical distortion included.  If False, this
            method will expect TAN_PIXEL coordinates, which are the pixel
            coordinates with estimated optical distortion removed. See the
            documentation in afw.cameraGeom for more details. (the default is
            True.)

        Returns
        -------
        numpy.ndarray
            A 2-D numpy array in which the first row is the RA coordinate and
            the second row is the Dec coordinate (both in degrees; in the
            International Celestial Reference System).
        """

        if isinstance(chipName, np.ndarray):
            chipNameList = chipName.tolist()
        else:
            chipNameList = chipName

        return raDecFromPixelCoords(
            xPix,
            yPix,
            chipNameList,
            camera=self._camera,
            obs_metadata=self._obs,
            epoch=epoch,
            includeDistortion=includeDistortion,
        )

    def pixelCoordsFromRaDec(
        self, ra, dec, chipName=None, epoch=2000.0, includeDistortion=True
    ):
        """Get the pixel positions (or nan if not on a chip) for objects based
        on their RA, and Dec (in degrees).

        Parameters
        ----------
        ra : float or numpy.ndarray
            ra is in degrees in the International Celestial Reference System.
        dec : float or numpy.ndarray
            dec is in degrees in the International Celestial Reference System.
        chipName : numpy.ndarray, str, None
            chipName designates the names of the chips on which the pixel
            coordinates will be reckoned. If an array, there must be as many
            chipNames as there are (RA, Dec) pairs. If a single value, all of
            the pixel coordinates will be reckoned on the same chip. If None,
            this method will calculate which chip each(RA, Dec) pair actually
            falls on, and return pixel coordinates for each (RA, Dec) pair on
            the appropriate chip. (the default is None.)
        epoch : float, optional
            epoch is the mean epoch in years of the celestial coordinate
            system. (the default is 2000.0.)
        includeDistortion : bool, optional
            If True (default), then this method will expect the true pixel
            coordinates with optical distortion included.  If False, this
            method will expect TAN_PIXEL coordinates, which are the pixel
            coordinates with estimated optical distortion removed. See the
            documentation in afw.cameraGeom for more details. (the default is
            True.)

        Returns
        -------
        numpy.ndarray
            A 2-D numpy array in which the first row is the x pixel coordinate
            and the second row is the y pixel coordinate.
        """

        return pixelCoordsFromRaDec(
            ra,
            dec,
            obs_metadata=self._obs,
            chipName=chipName,
            camera=self._camera,
            epoch=epoch,
            includeDistortion=includeDistortion,
        )

    def focalPlaneCoordsFromRaDec(self, ra, dec, epoch=2000.0):
        """Get the focal plane coordinates for all objects in the catalog.

        Parameters
        ----------
        ra : float or numpy.ndarray
            ra is in degrees in the International Celestial Reference System.
        dec : float or numpy.ndarray
            dec is in degrees in the International Celestial Reference System.
        epoch : float, optional
            epoch is the mean epoch in years of the celestial coordinate
            system. (the default is 2000.0.)

        Returns
        -------
        numpy.ndarray
            A 2-D numpy array in which the first row is the x focal plane
            coordinate and the second row is the y focal plane coordinate
            (both in millimeters).
        """

        return focalPlaneCoordsFromRaDec(
            ra, dec, obs_metadata=self._obs, epoch=epoch, camera=self._camera
        )


if __name__ == "__main__":
    pass
