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
import pickle
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from numpy.ma import masked_array
from src.ioUtilsNoStitch import caseData
import matplotlib.patches as patches
from matplotlib import colors

from piv.src.ioUtilsNoStitch import caseData
from piv.src.miscTools import curl2D
from piv.src.miscTools import curl

with open(dataDir + "/pickle/p70Case6.pickle","rb") as handle:
    
    p70Case7 = pickle.load(handle)
                               
vort = p70Case7["phase0"]["Vortz"]
curlVort = curl2D(p70Case7["phase0"]["X"],p70Case7["phase0"]["Y"],p70Case7["phase0"]["Vx"],p70Case7["phase0"]["Vr"])

def computeFlux(x1, y1, XX, YY, velX, velY, qty):
    
    #x1, y1 - array of x, y coordinate values
    
    xLine1 = np.zeros(len(y1)) + x1[0]
    xLine2 = x1
    xLine3 = np.zeros(len(y1)) + x1[-1]
    xLine4 = np.flip(x1)
    
    xLine = np.concatenate((xLine1, xLine2, xLine3, xLine4))
    
    yLine1 = y1
    yLine2 = np.zeros(len(x1)) + y1[-1]
    yLine3 = np.flip(y1)
    yLine4 = np.zeros(len(x1)) + y1[0]
    
    yLine = np.concatenate((yLine1, yLine2, yLine3, yLine4))
    
    #Midpoint for flux calculation
    line1 = np.stack((0.5*(xLine1[:-1] + xLine1[1:]),0.5*(yLine1[:-1] + yLine1[1:])),axis=-1)
    line2 = np.stack((0.5*(xLine2[:-1] + xLine2[1:]),0.5*(yLine2[:-1] + yLine2[1:])),axis=-1)
    line3 = np.stack((0.5*(xLine3[:-1] + xLine3[1:]),0.5*(yLine3[:-1] + yLine3[1:])),axis=-1)
    line4 = np.stack((0.5*(xLine4[:-1] + xLine4[1:]),0.5*(yLine4[:-1] + yLine4[1:])),axis=-1)
    line = np.concatenate((line1,line2,line3,line4),axis=0)
    
    X = XX[0,:]
    Y = YY[:,0]
    
    velXIntp = RegularGridInterpolator((X,Y),np.transpose(velX),method="cubic")
    velYIntp = RegularGridInterpolator((X,Y),np.transpose(velY),method="cubic")
    qtyIntp = RegularGridInterpolator((X,Y),np.transpose(qty),method="cubic")
    
    qtyFluxLine1 = qtyIntp(line1) * velXIntp(line1)
    qtyFluxLine2 = qtyIntp(line2) * velYIntp(line2)
    qtyFluxLine3 = qtyIntp(line3) * velXIntp(line3)
    qtyFluxLine4 = qtyIntp(line4) * velYIntp(line4)
    
    qtyFlowLine1 = qtyFluxLine1 * -np.diff(yLine1)#-ve symbol for left facing normal vector to be negative and right facing to be positive
    qtyFlowLine2 = qtyFluxLine2 * np.diff(xLine2)
    qtyFlowLine3 = qtyFluxLine3 * -np.diff(yLine3)
    qtyFlowLine4 = qtyFluxLine4 * np.diff(xLine4)
    
    massFlowLine1 = velXIntp(line1) * -np.diff(yLine1)
    massFlowLine2 = velYIntp(line2) * np.diff(xLine2)
    massFlowLine3 = velXIntp(line3) * -np.diff(yLine3)
    massFlowLine4 = velYIntp(line4) * np.diff(xLine4)
    
    qtyFlowLine = np.concatenate((qtyFlowLine1,qtyFlowLine2,qtyFlowLine3,qtyFlowLine4))
    massFlowLine = np.concatenate((massFlowLine1,massFlowLine2,massFlowLine3,massFlowLine4))
    
    print("Total vorticity flow:", np.sum(qtyFlowLine))
    print("Line1 vorticity flow:", np.sum(qtyFlowLine1))
    print("Line2 vorticity flow:", np.sum(qtyFlowLine2))
    print("Line3 vorticity flow:", np.sum(qtyFlowLine3))
    print("Line4 vorticity flow:", np.sum(qtyFlowLine4))
    print("\n")
    
    return xLine, yLine, line, qtyFlowLine, massFlowLine, {"line1" : qtyFlowLine1, "line2" : qtyFlowLine2, "line3" : qtyFlowLine3, "line4" : qtyFlowLine4}

XX = p70Case7["phase0"]["X"]
YY = p70Case7["phase0"]["Y"]
gridVelX = p70Case7["phase0"]["Vx"]
gridVelY = p70Case7["phase0"]["Vr"]
gridVortZ = p70Case7["phase0"]["Vortz"]
gridIsValid = p70Case7["phase0"]["ValidCells"]
maskedVelX = masked_array(gridVelX, mask=~gridIsValid)/-p70Case7["phase0"]["Vinf"]
maskedVortZ = masked_array(gridVortZ,mask=~gridIsValid)

yBounds = np.arange(0.07,0.16,0.001)

xLine, yLine, line, qtyFlowLine, massFlowLine, qtyFlowLineSplit= computeFlux(np.arange(-0.33,-0.30,0.001), yBounds, p70Case7["phase0"]["X"], p70Case7["phase0"]["Y"], p70Case7["phase0"]["Vx"],p70Case7["phase0"]["Vr"],p70Case7["phase0"]["Vortz"])
xLine, yLine, line, qtyFlowLine, massFlowLine, qtyFlowLineSplit= computeFlux(np.arange(-0.30,-0.27,0.001), yBounds, p70Case7["phase0"]["X"], p70Case7["phase0"]["Y"], p70Case7["phase0"]["Vx"],p70Case7["phase0"]["Vr"],p70Case7["phase0"]["Vortz"])
xLine, yLine, line, qtyFlowLine, massFlowLine, qtyFlowLineSplit= computeFlux(np.arange(-0.27,-0.24,0.001), yBounds, p70Case7["phase0"]["X"], p70Case7["phase0"]["Y"], p70Case7["phase0"]["Vx"],p70Case7["phase0"]["Vr"],p70Case7["phase0"]["Vortz"])
xLine, yLine, line, qtyFlowLine, massFlowLine, qtyFlowLineSplit= computeFlux(np.arange(-0.24,-0.21,0.001), yBounds, p70Case7["phase0"]["X"], p70Case7["phase0"]["Y"], p70Case7["phase0"]["Vx"],p70Case7["phase0"]["Vr"],p70Case7["phase0"]["Vortz"])
xLine, yLine, line, qtyFlowLine, massFlowLine, qtyFlowLineSplit= computeFlux(np.arange(-0.21,-0.18,0.001), yBounds, p70Case7["phase0"]["X"], p70Case7["phase0"]["Y"], p70Case7["phase0"]["Vx"],p70Case7["phase0"]["Vr"],p70Case7["phase0"]["Vortz"])
xLine, yLine, line, qtyFlowLine, massFlowLine, qtyFlowLineSplit= computeFlux(np.arange(-0.18,-0.15,0.001), yBounds, p70Case7["phase0"]["X"], p70Case7["phase0"]["Y"], p70Case7["phase0"]["Vx"],p70Case7["phase0"]["Vr"],p70Case7["phase0"]["Vortz"])
xLine, yLine, line, qtyFlowLine, massFlowLine, qtyFlowLineSplit= computeFlux(np.arange(-0.15,-0.12,0.001), yBounds, p70Case7["phase0"]["X"], p70Case7["phase0"]["Y"], p70Case7["phase0"]["Vx"],p70Case7["phase0"]["Vr"],p70Case7["phase0"]["Vortz"])
xLine, yLine, line, qtyFlowLine, massFlowLine, qtyFlowLineSplit= computeFlux(np.arange(-0.12,-0.09,0.001), yBounds, p70Case7["phase0"]["X"], p70Case7["phase0"]["Y"], p70Case7["phase0"]["Vx"],p70Case7["phase0"]["Vr"],p70Case7["phase0"]["Vortz"])
xLine, yLine, line, qtyFlowLine, massFlowLine, qtyFlowLineSplit= computeFlux(np.arange(-0.09,-0.06,0.001), yBounds, p70Case7["phase0"]["X"], p70Case7["phase0"]["Y"], p70Case7["phase0"]["Vx"],p70Case7["phase0"]["Vr"],p70Case7["phase0"]["Vortz"])

'''print("Total vorticity flow:", np.sum(qtyFlowLine))
print("Line1 vorticity flow:", np.sum(qtyFlowLineSplit["line1"]))
print("Line2 vorticity flow:", np.sum(qtyFlowLineSplit["line2"]))
print("Line3 vorticity flow:", np.sum(qtyFlowLineSplit["line3"]))
print("Line4 vorticity flow:", np.sum(qtyFlowLineSplit["line4"]))'''

divnorm = colors.TwoSlopeNorm(vmin=np.min(maskedVelX), vcenter=1.0, vmax=np.max(maskedVelX))

fig1, ax1 = plt.subplots()
#ax1.contourf(XX, YY, maskedVelX, cmap="seismic", norm=divnorm)
ax1.contourf(XX, YY, maskedVortZ, levels=np.linspace(-125,125,50), cmap="seismic")
ax1.plot(xLine, yLine, "g")
#ax1.axis("equal")
ax1.set_aspect("equal")
fig1.tight_layout(pad=0)


xmin, xmax = ax1.get_xlim()
ymin, ymax = ax1.get_ylim()
xy = (xmin,ymin)
width = xmax - xmin
height = ymax - ymin

#Patches
hatchedPatchVel = patches.Rectangle(xy, width, height, fc = "dimgrey", zorder=-10)
hatchedPatchVort = patches.Rectangle(xy, width, height, fc = "dimgrey", zorder=-10)
ax1.add_patch(hatchedPatchVel)

    
    
    
    