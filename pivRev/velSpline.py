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
from scipy.interpolate import CubicSpline
from scipy import signal as sig
from scipy.ndimage import median_filter

from src.ioUtilsPickle import pickleData

from src.miscTools import genPolygonLineBndry
from src.miscTools import genVector
from src.miscTools import nrmlVector2D
from src.miscTools import lowPassFilt

discPor = "45"
caseID = "00"

caseObj = pickleData(dataDir+"/pickledFilledData",discPor,caseID)
testVelExtp = caseObj.getVxExtp()
testCircRev = caseObj.getCirc2DRev()
testDiscVelRev = caseObj.getDiscVelRev()

circLine = caseObj.getCirc1D([-1.2,0],[1.7,0])[0]

#%%
plotPhaseID = "00"
phaseObj = getattr(caseObj,f"phase{plotPhaseID}")
maskedZ = np.ma.masked_array(phaseObj.gridWz,phaseObj.discMask)
maskedZv = np.ma.masked_array(phaseObj.gridVx,phaseObj.discMask)
maskedMask = np.ma.masked_array(phaseObj.discMask,np.logical_not(phaseObj.discMask))

discXc = phaseObj.discXc
discYc = phaseObj.discYc

# %%
# print("Circulation:", caseObj.circ1D[int(plotPhaseID)])
#Axial velocity plot
# x2 = circLine[int(plotPhaseID)][1,0]
# x1 = circLine[int(plotPhaseID)][0,0]

x1 = caseObj.gridX[0,15]
x2 = caseObj.gridX[0,-30]
y2 = circLine[int(plotPhaseID)][1,1]
y1 = circLine[int(plotPhaseID)][0,1]

lineDist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
minDist = caseObj.gridX[0,1] - caseObj.gridX[0,0]
nPts = int(lineDist/minDist) + 1
linePts = np.zeros((nPts,2))

linePts[:,0] = np.linspace(x1,x2,nPts)
linePts[:,1] = np.linspace(y1,y2,nPts)

VxLine = phaseObj.gridVxIntpr(linePts)
VyLine =  phaseObj.gridVyIntpr(linePts)
VxFftFilt = lowPassFilt(VxLine, len(VxLine), 20)

VxMinIndxs = sig.argrelmin(VxLine)[0]
VxMaxIndxs = sig.argrelmax(VxLine)[0]

# x3 = x1 + (x2 - x1)*0.3167
# x4 = x1 + (x2 - x1)*0.567

x3 = discXc - 0.25
x4 = discXc + 0.5

print("x3 wrt xc:", discXc - x3)
print("x4 wrt xc:", x4 - discXc)

x3Indx = np.argmin(np.abs(linePts[:,0] - x3))
x4Indx = np.argmin(np.abs(linePts[:,0] - x4))

line1Dist = np.sqrt((x3 - x1)**2 + (y2 - y1)**2)
line2Dist = np.sqrt((x2 - x4)**2 + (y2 - y1)**2)

n1Pts = int(line1Dist/minDist) + 1
n2Pts = int(line2Dist/minDist) + 1

line1Pts = np.zeros((n1Pts,2))
line2Pts = np.zeros((n2Pts,2))

line1Pts[:,0] = np.linspace(x1,x3,n1Pts)
line1Pts[:,1] = np.linspace(y1,y2,n1Pts)
line2Pts[:,0] = np.linspace(x4,x2,n2Pts)
line2Pts[:,1] = np.linspace(y1,y2,n2Pts)

line12Pts = np.row_stack((line1Pts,line2Pts))
VxLine12 = phaseObj.gridVxIntpr(line12Pts)
VyLine12 =  phaseObj.gridVyIntpr(line12Pts)
VxFftFilt12 = np.concatenate((VxFftFilt[:x3Indx+1],VxFftFilt[x4Indx:]))

VxSpline12 = CubicSpline(line12Pts[:,0], VxLine12)(linePts[:,0])
VxFiltSpline12 = CubicSpline(np.concatenate((linePts[:x3Indx+1,0],linePts[x4Indx:,0])), VxFftFilt12)(linePts[:,0])
VxFilt = median_filter(VxLine,100)

fig1, ax1 = plt.subplots()
ax1.plot(linePts[:,0] - discXc, VxLine, label=r"$V_x^{*}\ raw$")
# ax1.plot(linePts[:,0] - discXc, VxSpline12)
# ax1.plot(linePts[:,0] - discXc, VxFiltSpline12)
ax1.plot(testVelExtp[0][:,0] - discXc, testVelExtp[1][int(plotPhaseID)](testVelExtp[0][:,0]), label=r"$V_x^{*}\ extp$")
# ax1.plot(linePts[:,0] - discXc, VxFftFilt)
# ax1.plot(linePts[:,0] - discXc, VxFilt)
ax1.axvline(line1Pts[-1,0] - discXc,color="k",linestyle="--")
ax1.axvline(line2Pts[0,0] - discXc,color="k",linestyle="--")
ax1.set_xlabel(r"$x^{*}\ [-]$")
ax1.set_ylabel(r"$V_x^{*}\ [-]$")
ax1.legend()


fig,ax = plt.subplots()
# ax.contourf(caseObj.gridX, caseObj.gridY, maskedZ, levels = np.linspace(-100,100,101), extend="both")
ax.contourf(caseObj.gridX - discXc, caseObj.gridY - discYc, maskedMask, cmap="gray")
velPlot = ax.contourf(caseObj.gridX - discXc, caseObj.gridY - discYc, maskedZv, levels = np.linspace(0,np.max(maskedZv),31))
ax.plot(linePts[:,0] - discXc,linePts[:,1] - discYc,"k")
ax.axvline(line1Pts[-1,0] - discXc,color="k",linestyle="--")
ax.axvline(line2Pts[0,0] - discXc,color="k",linestyle="--")
ax.set_aspect("equal")
ax.set_xlabel(r"$x^{*}\ [-]$")
ax.set_ylabel(r"$y^{*}\ [-]$")
velCbar = fig.colorbar(velPlot, orientation="horizontal",aspect=40)
velCbar.set_label(r"$V_x^{*}\ [-]$")