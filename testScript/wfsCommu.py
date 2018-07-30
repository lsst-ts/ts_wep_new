import time
import numpy as np

from lsst.ts.wep.WEPController import WEPController

if __name__ == "__main__":

    # This script is to show the communication of WEP controller by SAL
    
    # Initiate the WEP Controller
    wepCntlr = WEPController()

    # Set the middle ware
    topicList = ["wavefrontErrorCalculated", "wavefrontError"]
    wepCntlr.setMiddleWare(topicList)

    # Issue the event
    sensorId = 11
    timestamp = time.time()
    priority = 1
    eventData = {"sensorId": sensorId,
                 "timestamp": timestamp,
                 "priority": priority}
    wepCntlr.issueEvent(topicList[0], eventData)
    time.sleep(1)

    # Issue the telemetry
    zkList = np.random.rand(19)
    telData = {"sensorId": sensorId,
               "annularZerikePolynomials": zkList,
               "timestamp": timestamp}
    wepCntlr.issueTelemetry(topicList[1], telData)
    time.sleep(10)

    # Turn off the SAL
    wepCntlr.shutDownSal()