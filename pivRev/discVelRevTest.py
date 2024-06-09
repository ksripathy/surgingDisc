#%%
import os
import sys
import pickle
import math

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

import numpy as np
import matplotlib.pyplot as plt
from src.ioUtilsPickle import pickleData

discPor = "70"
caseID = "05"
phaseID = "00"
yOffset = 0

caseObj = pickleData(dataDir+"/pickledFilledNdimData",discPor,caseID)
phaseObj = getattr(caseObj,f"phase{phaseID}")
VxExtpObj =caseObj.getVxExtp()[1]
avgVelTest = caseObj.getAvgDiscVel()
avgVelTest2 = caseObj.getAvgDiscVelV2(offset=0.2,annuliQty=10)
discVelTest = caseObj.getDiscVelRev()
polyObjs, discVelTestOld = caseObj.getDiscVelExtpV2(offset=0.2)

fig2,ax2 = plt.subplots()
# ax2.plot(avgVelTest)
ax2.plot(avgVelTest2,marker="x")
# ax2.plot(discVelTest[10,:])
ax2.plot(discVelTestOld[10,:],marker="x")
# %% Convergence analysis

cutOffs2 = np.arange(0.25,0.76,0.01)
VxExtp = np.zeros((len(VxExtpObj),len(cutOffs2)))

for i in range(len(VxExtp)):
    
    discXc = getattr(caseObj,caseObj.phaseObjNames[i]).discXc
    
    for j in range(len(cutOffs2)):
        
        VxExtpObj =caseObj.getVxExtp(yOffset=yOffset,cutOff2=cutOffs2[j])[1]
        if i==int(phaseID) and j == 25:
            plotExtpObj = VxExtpObj[i]
        VxExtp[i,j] = VxExtpObj[i](discXc)
        
xLocs = caseObj.gridX[0,30:-30]
yLoc = caseObj.phase00.discYc + yOffset
pts = np.column_stack((xLocs,yLoc*np.ones(len(xLocs))))

xEndIndex = np.argmin(np.abs(caseObj.gridX[0,:] - phaseObj.discXc))
xLocs2 = caseObj.gridX[0,30:xEndIndex+1]

xPlot = xLocs - phaseObj.discXc
xPlot2 = xLocs2 - phaseObj.discXc
VxRawPlot = phaseObj.gridVxIntpr(pts)
VxExtpPlot = plotExtpObj(xLocs)

fig1, ax1 = plt.subplots()
ax1.plot(xPlot,VxRawPlot)
ax1.plot(xPlot,VxExtpPlot)
ax1.plot(xPlot2,polyObjs[10,int(phaseID)](xLocs2))



