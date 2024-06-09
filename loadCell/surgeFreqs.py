import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from src.miscTools import harmonicPassFilt
from src.miscTools import freqsSignal
from src.miscTools import cycleAvgStd

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

from src.recordedPhasesForcesRev import recordedData

discPor = "45"
loadCell = "20N"
caseNo="05"

logFilePath = dataDir + f"/testMatrixTotal/curved/{loadCell}/p{discPor}LogCase{caseNo}.txt"
totalLoadsFilePath = dataDir + f"/testMatrixTotal/curved/{loadCell}/p{discPor}LoadCase{caseNo}.txt"

caseObj = recordedData(totalLoadsFilePath,logFilePath)
freq = caseObj.surgeFreq