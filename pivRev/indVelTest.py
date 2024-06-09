#%%

import os
import sys
import numpy as np
import matplotlib.pyplot as plt


#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

from src.pivFields import planarPIVField
from src.miscTools import curl2D
from src.miscTools import bool2DMinMax

p45DiscAxialLocs = np.array([[3.855, -19.480, -27.204, -17.344, 7.470, 40.008, 66.630, 75.350, 63.672, 35.407, 3.362, -19.973],
                               [12.236, 0.733, -2.554, 4.512, 18.316, 35.078, 47.403, 51.018, 43.130, 28.669, 11.743, 0.897],
                               [1.883, -21.617, -27.368, -18.659, 7.470, 40.665, 70.081, 75.997, 65.808, 35.900, 1.390, -22.110],
                               [-3.075, -36.559, -49.008, -38.735, -3.568, 45.403, 89.116, 103.740, 88.294, 46.882, -4.390, -38.571],
                               [3.082, -20.322, -25.685, -16.909, 8.932, 41.600, 71.342, 77.192, 66.953, 36.724, 2.431, -20.322],
                               [-3.897, -42.022, -47.166, -39.064, -4.390, 46.553, 93.060, 107.190, 90.266, 46.882, -5.047, -41.693],
                               [25.683, None, None, None, None, None, None, None, None, None, None, None],
                               [25.683, None, None, None, None, None, None, None, None, None, None, None]])

p70DiscAxialLocs = np.array([[-0.742, -24.672, -31.884, -22.213, 2.536, 34.498, 61.214, 69.901, 57.936, 29.908, -1.397, -24.836],
                               [6.962, -4.184, -7.134, -0.414, 13.026, 30.072, 41.710, 45.807, 37.776, 23.680, 6.634, -4.012],
                               [-3.364, -27.131, -31.720, -23.852, 2.372, 35.153, 64.165, 69.246, 60.559, 29.581, -3.692, -27.294],
                               [-5.207, -38.687, -52.827, -38.687, -4.395, 43.225, 88.082, 103.20, 86.781, 43.225, -6.020, -39.337],
                               [3.082, -20.322, -25.685, -16.909, 8.932, 41.600, 71.342, 77.192, 66.953, 36.724, 2.431, -20.322],
                               [-6.680, -43.328, -50.876, -40.312, -5.857, 44.850, 90.682, 105.150, 87.269, 44.688, -7.645, -42.263],
                               [26.012, None, None, None, None, None, None, None, None, None, None, None],
                               [26.012, None, None, None, None, None, None, None, None, None, None, None]])

p45DiscCentreLocs = np.array([5.932, 5.932, 5.932, 8.069, 8.069, 8.069, 8.069, 8.069])
p70DiscCentreLocs = np.array([8.069, 8.069, 8.069, 6.369, 6.369, 6.369, 8.069, 8.069])

discPor = "45"

if discPor == "45":
    
    discAxialLocs = p45DiscAxialLocs
    discCentreLocs = p45DiscCentreLocs
    
elif discPor == "70":
    
    discAxialLocs = p70DiscAxialLocs
    discCentreLocs = p70DiscCentreLocs

caseID = "06"
caseNo = int(caseID)

phaseID = "00"
phaseNo = int(phaseID)

discXc = discAxialLocs[caseNo,phaseNo]
discYc = discCentreLocs[caseNo]

fieldObj = planarPIVField(dataDir + f"/p{discPor}MeanCase{caseID}Phase{phaseID}.dat", dataDir + f"/p{discPor}LogCase{caseID}Phase{phaseID}.txt", discXc, discYc)

ustrmOffsets = np.linspace(0.05,0.5,10)

fig, ax = plt.subplots()

fieldObj.getIndVelRaw(fieldObj.discUstreamLine)
indVxu = fieldObj.indVxRaw
fieldObj.getIndVelRaw(fieldObj.discCentreLine)
indVxc = fieldObj.indVxRaw

ax.plot(fieldObj.discUstreamLine[:,1],indVxu)
ax.plot(fieldObj.discUstreamLine[:,1],indVxc)

#%%

lineTemplate = np.column_stack((np.ones(21),np.zeros(21)))*fieldObj.ndimField.discDia

fieldObj.getIndVelIntp(fieldObj.discUstreamLine)
indVxu1RBF = fieldObj.indVxIntp

tempLine2 = 0.05 * lineTemplate
fieldObj.getIndVelIntp(fieldObj.discUstreamLine + tempLine2)
indVxu2RBF = fieldObj.indVxIntp

fieldObj.getIndVelIntp(fieldObj.discCentreLine)
indVxcRBF = fieldObj.indVxIntp

tempLine3 = -0.05 * lineTemplate
fieldObj.getIndVelIntp(fieldObj.discUstreamLine + tempLine3)
indVxu3RBF = fieldObj.indVxIntp

tempLine4 = -0.10 * lineTemplate
fieldObj.getIndVelIntp(fieldObj.discUstreamLine + tempLine4)
indVxu4RBF = fieldObj.indVxIntp

tempLine5 = -0.15 * lineTemplate
fieldObj.getIndVelIntp(fieldObj.discUstreamLine + tempLine5)
indVxu5RBF = fieldObj.indVxIntp

tempLine6 = -0.40 * lineTemplate
fieldObj.getIndVelIntp(fieldObj.discUstreamLine + tempLine6)
indVxu6RBF = fieldObj.indVxIntp

fig2, ax2 = plt.subplots()
ax2.plot(fieldObj.discUstreamLine[:,1],indVxcRBF)
ax2.plot(fieldObj.discUstreamLine[:,1],indVxu2RBF)
ax2.plot(fieldObj.discUstreamLine[:,1],indVxu1RBF)
ax2.plot(fieldObj.discUstreamLine[:,1],indVxu,"k--")
ax2.plot(fieldObj.discUstreamLine[:,1],indVxu3RBF)
ax2.plot(fieldObj.discUstreamLine[:,1],indVxu4RBF)
ax2.plot(fieldObj.discUstreamLine[:,1],indVxu5RBF)
ax2.plot(fieldObj.discUstreamLine[:,1],indVxu6RBF)


# for offset in ustrmOffsets:
    
#     fieldObj.getIndVelRaw(ustrOffset=offset)
#     plotx = fieldObj.lineRaw[:,1]
#     ploty = fieldObj.indVxRaw
#     ax.plot(plotx, ploty, label=f"{offset}")
# %%
