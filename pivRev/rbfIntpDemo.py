#%%
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob


#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

from src.ioUtils import caseData
from src.pivFields import planarPIVField
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
caseID = "03"
phaseID = "11"

if discPor == "45":
    
    discAxialLocs = p45DiscAxialLocs
    discCentreLocs = p45DiscCentreLocs
    
elif discPor == "70":
            
    discAxialLocs = p70DiscAxialLocs
    discCentreLocs = p70DiscCentreLocs
            
discXc = discAxialLocs[int(caseID),int(phaseID)]
discYc = discCentreLocs[int(caseID)]

caseDataPath = sorted(glob.glob(dataDir + f"/p{discPor}MeanCase{caseID}Phase{phaseID}*"))[0]
caseLogPath = sorted(glob.glob(dataDir + f"/p{discPor}LogCase{caseID}Phase{phaseID}*"))[0]

# caseObj = caseData(dataDir,"45","06")
phaseObj = planarPIVField(caseDataPath, caseLogPath, discXc,discYc)
maskedField = phaseObj.frame1.maskFrame()

#%%
# rbfMask = phaseObj.getRBFMask(phaseObj.ndimField)
# minY, maxY, minX, maxX = bool2DMinMax(rbfMask)

rbfMask1, rbfMask2 = phaseObj.getRBFMask3(phaseObj.ndimField,0.1,0.16,1,1)[:2]

ustrMinY, ustrMaxY, ustrMinX, ustrMaxX = bool2DMinMax(rbfMask1)
dstrMinY, dstrMaxY, dstrMinX, dstrMaxX = bool2DMinMax(rbfMask2)

fig, ax = plt.subplots()
ax.contourf(maskedField.gridX, maskedField.gridY, maskedField.gridWz,cmap="seismic",levels=np.linspace(-75,75,150),extend="both")
ax.set_aspect("equal")

# mask = np.array(maskedField.gridX[minY:maxY,minX:maxX])
# mask[:,:] = 100

mask1 = np.array(maskedField.gridX[ustrMinY:ustrMaxY,ustrMinX:ustrMaxX])
mask1[:,:] = 100
mask2 = np.array(maskedField.gridX[dstrMinY:dstrMaxY,dstrMinX:dstrMaxX])
mask2[:,:] = 100

fig2, ax2 = plt.subplots()
ax2.contourf(maskedField.gridX, maskedField.gridY, maskedField.gridWz,cmap="seismic",levels=np.linspace(-75,75,150),extend="both")
# ax2.contourf(maskedField.gridX[minY:maxY,minX:maxX], maskedField.gridY[minY:maxY,minX:maxX],mask)
ax2.contourf(maskedField.gridX[ustrMinY:ustrMaxY,ustrMinX:ustrMaxX], maskedField.gridY[ustrMinY:ustrMaxY,ustrMinX:ustrMaxX],mask1)
ax2.contourf(maskedField.gridX[dstrMinY:dstrMaxY,dstrMinX:dstrMaxX], maskedField.gridY[dstrMinY:dstrMaxY,dstrMinX:dstrMaxX],mask2)
ax2.set_aspect("equal")


# %%
