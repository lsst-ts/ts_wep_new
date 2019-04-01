import os
import numpy as np

from lsst.ts.wep.cwfs.Tool import getConfigValue


class Instrument(object):

    def __init__(self, instruFolder):
        """Instrument class for wavefront estimation.

        Parameters
        ----------
        instruFolder : str
            Path to instrument folder.
        """

        # Folder of instrument
        self.instDir = instruFolder

        # Instrument parameters
        self.filename = None
        self.maskParam = None
        self.parameter = None
        self.xSensor = None
        self.ySensor = None
        self.xoSensor = None
        self.yoSensor = None
        self.instName = None

    def getInstFileName(self):
        """Get the instrument file name.

        Returns
        -------
        str
            Instrument file name.
        """

        return self.filename

    def config(self, instruName, sensorSamples):
        """Configure the instrument.

        Parameters
        ----------
        instruName : str
            Instrument name. It is "lsst" in the baseline.
        sensorSamples : int
            Dimension of image on sensor.
        """

        # Name of instrument
        self.instName = instruName

        # Path of instrument file
        self.filename = os.path.join(self.instDir, instruName,
                                     (instruName + ".param"))

        # Directory to mask off-axis correction
        self.maskParam = os.path.join(self.instDir, instruName,
                                      "mask_migrate.txt")

        # Get the instrument parameters
        sensorSamples = int(sensorSamples)
        parameter = self._readFile(sensorSamples)
        self.parameter = parameter

        # Construct sensor coordinates
        sensorFactor = parameter["sensorFactor"]
        ySensor, xSensor = np.mgrid[-(sensorSamples/2-0.5):(sensorSamples/2 + 0.5),
                                    -(sensorSamples/2-0.5):(sensorSamples/2 + 0.5)]
        self.xSensor = xSensor/(sensorSamples/2/sensorFactor)
        self.ySensor = ySensor/(sensorSamples/2/sensorFactor)

        # Get the position index that is out of annular aperature range
        obscuration = parameter["obscuration"]
        r2Sensor = self.xSensor**2 + self.ySensor**2
        idx = (r2Sensor > 1) | (r2Sensor < obscuration**2)

        # o indicates annulus
        self.xoSensor = self.xSensor.copy()
        self.yoSensor = self.ySensor.copy()

        # Define the value to be NaN if it is not in pupul
        self.xoSensor[idx] = np.nan
        self.yoSensor[idx] = np.nan

    def _readFile(self, sensorSamples):
        """Read the instrument parameter file.

        Parameters
        ----------
        sensorSamples : int
            Dimension of image on sensor.

        Returns
        -------
        dict
            Instrument parameters.
        """

        # Get all parameters
        obscuration = getConfigValue(self.filename, "Obscuration", index=-1)
        focalLength = getConfigValue(self.filename, "Focal_length", index=-1)
        apertureDiameter = getConfigValue(self.filename, "Aperture_diameter",
                                          index=-1)
        offset = getConfigValue(self.filename, "Offset", index=-1)
        pixelSize = getConfigValue(self.filename, "Pixel_size", index=-1)

        # Marginal focal length (m = f/sqrt(f^2 - (D/2)^2))
        marginalFL = np.sqrt(focalLength**2 - (apertureDiameter/2)**2)

        # The below need to be instrument parameters, b/c it is not specific
        # for I1 or I2
        sensorFactor = sensorSamples/(offset*apertureDiameter/focalLength/pixelSize)

        # Collect all parameters of instrument into a single dictionary
        # attribute
        parameter = {"obscuration": obscuration,
                     "focalLength": focalLength,
                     "apertureDiameter": apertureDiameter,
                     "offset": offset,
                     "pixelSize": pixelSize,
                     "marginalFL": marginalFL,
                     "sensorSamples": sensorSamples,
                     "sensorFactor": sensorFactor}

        return parameter


if __name__ == "__main__":
    pass
