import time
import numpy as np

from lsst.ts.wep.WEPController import WEPController

if __name__ == "__main__":

    # This script is to show the communication of WEP controller by SAL
    
    # Initiate the WEP Controller
    wepCntlr = WEPController()

    # Set the middle ware
    topicList = ["WavefrontErrorCalculated", "WavefrontError"]
    wepCntlr.setMiddleWare(topicList)

    # Issue the event
    sensorName = "R:2,2 S:1,1"
    timestamp = time.time()
    priority = 1
    eventData = {"sensorID": sensorName,
                 "timestamp": timestamp,
                 "priority": priority}
    wepCntlr.issueEvent("WavefrontErrorCalculated", eventData)
    time.sleep(1)

    # Issue the telemetry
    zkList = np.random.rand(19)
    telData = {"sensorID": sensorName,
               "annularZerikePolynomials": zkList,
               "timestamp": timestamp}
    wepCntlr.issueTelemetry("WavefrontError", telData)
    time.sleep(10)

    # Turn off the SAL
    wepCntlr.shutDownSal()