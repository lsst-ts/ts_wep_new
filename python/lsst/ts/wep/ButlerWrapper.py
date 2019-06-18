import numpy as np
from lsst.daf.persistence import Butler


class ButlerWrapper(object):

    def __init__(self, inputs, outputs=None):
        """Initialize the butler wrapper class.

        Parameters
        ----------
        inputs : RepositoryArgs, dict, or str
            Can be a single item or a list. Provides arguments to load an
            existing repository (or repositories). String is assumed to be a
            URI and is used as the cfgRoot (URI to the location of the cfg
            file). (Local file system URI does not have to start with
            'file://' and in this way can be a relative path). The
            'RepositoryArgs' class can be used to provide more parameters with
            which to initialize a repository (such as 'mapper', 'mapperArgs',
            'tags', etc. See the 'RepositoryArgs' documentation for more
            details). A dict may be used as shorthand for a 'RepositoryArgs'
            class instance. The dict keys must match parameters to the
            'RepositoryArgs.__init__' function.
        outputs : RepositoryArgs, dict, or str, optional
            Provides arguments to load one or more existing repositories or
            create new ones. The different types are handled the same as for
            'inputs'. (the default is None.)
        """

        self._butler = Butler(inputs=inputs, outputs=outputs)

    def setInputsAndOutputs(self, inputs, outputs=None):
        """Set the inputs and outputs of butler.

        Parameters
        ----------
        inputs : RepositoryArgs, dict, or str
            Can be a single item or a list. Provides arguments to load an
            existing repository (or repositories). String is assumed to be a
            URI and is used as the cfgRoot (URI to the location of the cfg
            file). (Local file system URI does not have to start with
            'file://' and in this way can be a relative path). The
            'RepositoryArgs' class can be used to provide more parameters with
            which to initialize a repository (such as 'mapper', 'mapperArgs',
            'tags', etc. See the 'RepositoryArgs' documentation for more
            details). A dict may be used as shorthand for a 'RepositoryArgs'
            class instance. The dict keys must match parameters to the
            'RepositoryArgs.__init__' function.
        outputs : RepositoryArgs, dict, or str, optional
            Provides arguments to load one or more existing repositories or
            create new ones. The different types are handled the same as for
            'inputs'. (the default is None.)
        """

        self._butler = Butler(inputs=inputs, outputs=outputs)

    def getRawExp(self, visit, raft, sensor, snap=None):
        """Get the raw exposure.

        Parameters
        ----------
        visit : int
            Visit Id.
        raft : str
            Abbreviated raft name (e.g. "R22").
        sensor : str
            Abbreviated sensor name (e.g. "S11").
        snap : int, optional
            Snap time (0 or 1) means first/ second exposure. (the default is
            None.)

        Returns
        -------
        lsst.afw.image.exposure.exposure.ExposureF
            Raw exposure object.
        """

        # Data Id keys in w_2019_24:
        # {'run': str, 'raftName': str, 'visit': int, 'detectorName': str,
        #  'detector': int, 'snap': int}

        dataId = self._getDefaultDataId(visit, raft, sensor)
        self._extendDataId(dataId, snap=snap)

        return self._butler.get("raw", dataId=dataId)

    def _getDefaultDataId(self, visit, raft, sensor):
        """Get the default data Id.

        Parameters
        ----------
        visit : int
            Visit Id.
        raft : str
            Abbreviated raft name (e.g. "R22").
        sensor : str
            Abbreviated sensor name (e.g. "S11").

        Returns
        -------
        dict
            Default data Id.
        """

        dataId = dict(visit=int(visit), raftName=raft, detectorName=sensor)

        return dataId

    def _extendDataId(self, dataId, snap=None, aFilter=None):
        """Extend the data Id.

        Parameters
        ----------
        dataId : dict
            Data Id.
        snap : int, optional
            Snap time (0 or 1) means first/ second exposure. (the default is
            None.)
        aFilter : str, optional
            Active filter ("u", "g", "r", "i", "z", "y") (the default is None.)
        """

        if (snap is not None) and isinstance(snap, (int, float)):
            dataId["snap"] = int(snap)

        if (aFilter is not None) and isinstance(aFilter, str):
            dataId["filter"] = aFilter

    def getPostIsrCcd(self, visit, raft, sensor, aFilter=None):
        """Get the post-ISR CCD exposure.

        ISR: Instrument signature removal.
        CCD: Charge-coupled device.

        Parameters
        ----------
        visit : int
            Visit Id.
        raft : str
            Abbreviated raft name (e.g. "R22").
        sensor : str
            Abbreviated sensor name (e.g. "S11").
        aFilter : str, optional
            Active filter ("u", "g", "r", "i", "z", "y") (the default is None.)

        Returns
        -------
        lsst.afw.image.exposure.exposure.ExposureF
            Post-ISR CCD object.
        """

        # Data Id keys in w_2019_24:
        # {'visit': int, 'filter': str, 'raftName': str, 'detectorName': str,
        #  'detector': int}

        dataId = self._getDefaultDataId(visit, raft, sensor)
        self._extendDataId(dataId, aFilter=aFilter)

        return self._butler.get("postISRCCD", dataId=dataId)

    def getEimage(self, visit, raft, sensor, snap=None):
        """Get the PhoSim eimage exposure.

        Parameters
        ----------
        visit : int
            Visit Id.
        raft : str
            Abbreviated raft name (e.g. "R22").
        sensor : str
            Abbreviated sensor name (e.g. "S11").
        snap : int, optional
            Snap time (0 or 1) means first/ second exposure. (the default is
            None.)

        Returns
        -------
        lsst.afw.image.exposure.exposure.ExposureF
            Eimage exposure object.
        """

        # Data Id keys in w_2019_24:
        # {'visit': int, 'snap': int, 'raftName': str, 'detectorName': str}
        dataId = self._getDefaultDataId(visit, raft, sensor)
        self._extendDataId(dataId, snap=snap)

        return self._butler.get("eimage", dataId=dataId)

    @staticmethod
    def getImageData(exposure):
        """Get the image data.

        Parameters
        ----------
        exposure : lsst.afw.image.exposure.exposure.ExposureF
            Exposure object.

        Returns
        -------
        numpy.ndarray
            Image data.
        """

        # Get the numpy array data based on the input object type
        if isinstance(exposure, np.ndarray):
            data = exposure
        elif hasattr(exposure, "getMaskedImage"):
            data = exposure.getMaskedImage().getImage().getArray()
        elif hasattr(exposure, "getImage"):
            data = exposure.getImage().getArray()
        else:
            data = exposure.getArray()

        # Return the data in numpy array
        return data


if __name__ == "__main__":
    pass
