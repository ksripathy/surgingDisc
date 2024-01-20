import os
import sys
#%matplotlib widget

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

#Add src folder to python path
sys.path.append("..")

import numpy as np
import matplotlib.pyplot as plt

from piv.src.ioUtilsNoStitch import caseData

p70Case0 = caseData(dataDir + "/tecPlotNoStitch/p70Case0*")

