import os
from pathlib import Path
import lsst.ts.wep

def getModulePath():
    """
    
    Get the directory of WEP module.
    
    Returns:
        [str] -- Directory of WEP module.
    """

    # Get the path of module
    modulePathList = os.path.dirname(lsst.ts.wep.__file__).split(os.sep)[3:-4]
    modulePath = os.path.join(str(Path.home()), *modulePathList)
    
    return modulePath

if __name__ == "__main__":
    pass