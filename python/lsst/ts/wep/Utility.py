import os
import lsst.ts.wep


def getModulePath(module=lsst.ts.wep, startIdx=1, endIdx=-4):
    """Get the path of module.
    Parameters
    ----------
    module : str, optional
        Module name. (the default is lsst.ts.ofc.)
    startIdx : int, optional
        Start index. (the default is 1.)
    endIdx : int, optional
        End index. (the default is -4.)
    Returns
    -------
    str
        Directory path of module based on the start and end indexes.
    """

    # Get the path of module
    modulePathList = os.path.dirname(module.__file__).split(
                                os.sep)[int(startIdx):int(endIdx)]
    modulePath = os.path.join(os.sep, *modulePathList)

    return modulePath


def runProgram(command, binDir=None, argstring=None):
    """Run the program w/o arguments.

    Parameters
    ----------
    command : str
        Command of application.
    binDir : str, optional
        Directory of binary application. (the default is None.)
    argstring : str, optional
        Arguments of program. (the default is None.)

    Raises
    ------
    RuntimeError
        Error running of command.
    """

    # Directory of binary application
    if (binDir is not None):
        command = os.path.join(binDir, command)

    # Arguments for the program
    if (argstring is not None):
        command += (" " + argstring)

    # Call the program w/o arguments
    if (subprocess.call(command, shell=True) != 0):
        raise RuntimeError("Error running: %s" % command)


if __name__ == "__main__":
    pass
