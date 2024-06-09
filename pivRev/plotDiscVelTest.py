import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

from src.plotDiscVel import plotDiscVel

plotObj = plotDiscVel(dataDir+"/pickledPlotData",plotDir+"/fieldsAeroDays")
# plotObj.pcolorPlot("70","05",discVelAttrb="discVxcExtp")
# plotObj.casesPlot("70",discVelAttrb="discVxcExtp",spanID="10")
# plotObj.staticPlot("70",discVelAttrb="discVxcExtp")
# plotObj.dynInflowPlot("70",discVelAttrb="discVxcExtp")
plotObj.fieldPlot("70","07")