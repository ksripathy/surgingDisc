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
    
# points = np.array([[0,0],[1,0],[1,1],[0,1]])
# points = np.array([[0,0],[1,2]])

# lineSegs, cntrlPtSegs = genPolygonLineBndry(points,0.1)
# vec1 = genVector(lineSegs[0][:-1,:],lineSegs[0][1:,:])
# nrmlVec1 = nrmlVector2D(lineSegs[0][:-1,:],lineSegs[0][1:,:])[1]

discPor = "70"
caseID = "05"

caseObj = pickleData(dataDir+"/pickledFilledData",discPor,caseID)

circLine = caseObj.getCirc1D([-0.25,0],[1,0],False)
circLine2 = caseObj.getCirc1D([-0.25,0],[1,0],True)
circLine3 = caseObj.getCirc2DRev()[1]

fig2, ax2 = plt.subplots()
ax2.plot(circLine[1])
ax2.plot(circLine2[1])
ax2.plot(-circLine3[:,2])

#%%

plotPhaseID = "07"
phaseObj = getattr(caseObj,f"phase{plotPhaseID}")
maskedZ = np.ma.masked_array(phaseObj.gridWz,phaseObj.discMask)
maskedZv = np.ma.masked_array(phaseObj.gridVx,phaseObj.discMask)
maskedMask = np.ma.masked_array(phaseObj.discMask,np.logical_not(phaseObj.discMask))

discXc = phaseObj.discXc
discYc = phaseObj.discYc

fig,ax = plt.subplots()
# ax.contourf(caseObj.gridX, caseObj.gridY, maskedZ, levels = np.linspace(-100,100,101), extend="both")
ax.contourf(caseObj.gridX, caseObj.gridY, maskedMask, cmap="gray")
ax.contourf(caseObj.gridX, caseObj.gridY, maskedZv, levels = np.linspace(0,np.max(maskedZv),31))
ax.plot(circLine[int(plotPhaseID)][:,0],circLine[int(plotPhaseID)][:,1],"k")
ax.set_aspect("equal")



# %%
print("Circulation:", caseObj.circ1D[int(plotPhaseID)])
#Axial velocity plot
x2 = circLine[int(plotPhaseID)][1,0]
x1 = circLine[int(plotPhaseID)][0,0]
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

x3 = x1 + (x2 - x1)*0.35
x4 = x1 + (x2 - x1)*0.60
# x4OldIndx = np.argmin(np.abs(linePts[:,0] - x4Old))

# VxMinIndxs = sig.argrelmin(VxLine,order=20)[0]
# VxMaxIndxs = sig.argrelmax(VxLine,order=50)[0]
# x4NewIndx = VxMinIndxs[np.argmin(np.abs(VxMinIndxs - x4OldIndx))]
# x4 = linePts[x4NewIndx,0]

line1Dist = np.sqrt((x3 - x1)**2 + (y2 - y1)**2)
line2Dist = np.sqrt((x2 - x4)**2 + (y2 - y1)**2)
line3Dist = np.sqrt((x4 - x3)**2 + (y2 - y1)**2)
n1Pts = int(line1Dist/minDist) + 1
n2Pts = int(line2Dist/minDist) + 1
n3Pts = int(line3Dist/minDist) + 1

line1Pts = np.zeros((n1Pts,2))
line2Pts = np.zeros((n2Pts,2))
line3Pts = np.zeros((n3Pts,2))

line1Pts[:,0] = np.linspace(x1,x3,n1Pts)
line1Pts[:,1] = np.linspace(y1,y2,n1Pts)
VxLine1 = phaseObj.gridVxIntpr(line1Pts)
VyLine1 = phaseObj.gridVyIntpr(line1Pts)

line2Pts[:,0] = np.linspace(x4,x2,n2Pts)
line2Pts[:,1] = np.linspace(y1,y2,n2Pts)

line3Pts[:,0] = np.linspace(x4,x3,n3Pts)
line3Pts[:,1] = np.linspace(y1,y2,n3Pts)

line12Pts = np.row_stack((line1Pts,line2Pts))
VxLine12 = phaseObj.gridVxIntpr(line12Pts)
VyLine12 =  phaseObj.gridVyIntpr(line12Pts)
VxLine3 = phaseObj.gridVxIntpr(line3Pts)


VxFit = np.poly1d(np.polyfit(linePts[:,0], VxLine, 3))(linePts[:,0])
VxFit1 = np.poly1d(np.polyfit(line1Pts[:,0], VxLine1, 3))(linePts[:,0])
VxFit12 = np.poly1d(np.polyfit(line12Pts[:,0], VxLine12, 3))(linePts[:,0])
VxSpline12 = CubicSpline(line12Pts[:,0], VxLine12)(linePts[:,0])
# VxSpline12v2 = make_interp_spline(line12Pts[:,0], VxLine12, k=5)(linePts[:,0])
VxFit3 = np.poly1d(np.polyfit(line3Pts[:,0], VxLine3, 3))(line3Pts[:,0])
VxFilt = median_filter(VxLine3,25,mode="nearest")
# VxFitFilt = np.poly1d(np.polyfit(linePts[:,0], VxFilt, 3))(linePts[:,0])

fig1, ax1 = plt.subplots()
ax1.plot(linePts[:,0] - discXc, VxLine)
# ax1.plot(linePts[:,0], VxFit)
# ax1.plot(linePts[:,0] - discXc, VxFit12)
# ax1.plot(line3Pts[:,0] - discXc, VxFit3)
ax1.plot(linePts[:,0] - discXc, VxSpline12)
# ax1.plot(line3Pts[:,0] - discXc, VxFilt)
# ax1.plot(linePts[:,0] - discXc, VxFitFilt)
# ax1.plot(linePts[:,0] - discXc, VxSpline12v2)
# ax1.plot(linePts[:,0] - discXc, VxFitSpline12)
# ax1.plot(linePts[:,0], VxFit1)
# ax1.plot(linePts[:,0], VxFit3, "--")
ax1.axvline(line1Pts[-1,0] - discXc,color="k",linestyle="--")
ax1.axvline(line2Pts[0,0] - discXc,color="k",linestyle="--")
# ax1.axvline(x4Old - discXc,color="r",linestyle="--")
# %%
