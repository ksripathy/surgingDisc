# %%
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

def stitchWeight(x, xMin, xMax, A1):
    
    xNorm = (x - xMin) / (xMax- xMin)
    
    return 1/(1 + np.exp(A1 * (xNorm - 0.5)))

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

#%%
x = fieldObj.ndimField.gridX
y = fieldObj.ndimField.gridY

maskOrigObj = fieldObj.ndimField.maskFrame()
maskFilledObj = fieldObj.ndimFieldFilled.maskFrame()

origZ = maskOrigObj.gridVx
filledZ = maskFilledObj.gridVx

figOrig, axOrig = plt.subplots()
axOrig.contourf(x,y,origZ)
axOrig.set_aspect("equal")

figFilled, axFilled = plt.subplots()
axFilled.contourf(x,y,filledZ)
axFilled.set_aspect("equal")


#%%
frameObj0 = fieldObj.frame0
maskFrameObj0 = frameObj0.maskFrame()
frameObj1 = fieldObj.frame1
maskFrameObj1 = frameObj1.maskFrame()

# frameOverlapIndxs = np.argwhere(np.logical_and(np.logical_not(frameObj0.gridMask),np.logical_not(frameObj1.gridMask)))
# # frameOVerlapIndxsFlattened = np.ravel_multi_index([frameOverlapIndxs[:0]])
# xMinIndx = min(frameOverlapIndxs[:,1])
# xMaxIndx = max(frameOverlapIndxs[:,1])

# xMin = frameObj0.gridX[0,xMinIndx]
# xMax = frameObj0.gridX[0,xMaxIndx]

# xNormLocs = (frameObj0.gridX[frameOverlapIndxs[:,0],frameOverlapIndxs[:,1]] - xMin)/(xMax - xMin)
# weights1 = stitchWeight(frameObj0.gridX[frameOverlapIndxs[:,0],frameOverlapIndxs[:,1]], xMin, xMax, 12)
# weights2 = stitchWeight(frameObj0.gridX[frameOverlapIndxs[:,0],frameOverlapIndxs[:,1]], xMin, xMax, -12)

# fig, ax = plt.subplots()
# ax.plot(xNormLocs[:91], weights1[:91])
# ax.plot(xNormLocs[:91], weights2[:91])

# shiftedFrame = frameObj.shiftFrame(frameObj, 0, -3)
# temp1 = frameObj.gridY
# temp2 = shiftedFrame.gridY

discLoc = discAxialLocs[caseNo,phaseNo]
flipDiscLoc = np.max(frameObj1.x) - (np.abs(np.min(frameObj1.x)) + discLoc)
discCentreLoc = discCentreLocs[caseNo]

combinedFrameObj = fieldObj.combineFrames(frameObj0,frameObj1,True)
flipCombinedFrameObj = combinedFrameObj.flipGrid()

maskCombinedFrameObj = combinedFrameObj.maskFrame()
maskFlipCombinedFrameObj = flipCombinedFrameObj.maskFrame()

# plotObj = maskCombinedFrameObj

# fig1, ax1 = plt.subplots()
# ax1.contourf(plotObj.gridX, plotObj.gridY, np.flip(plotObj.gridVx/(-3),axis=1), levels=np.linspace(-0.25,1.5,100))
# ax1.quiver(plotObj.gridX[::24,::24], plotObj.gridY[::24,::24],np.flip(-plotObj.gridVx,axis=1)[::24,::24],np.flip(plotObj.gridVy,axis=1)[::24,::24])
# ax1.plot(flipDiscLoc,3.959,"kx")
# ax1.set_aspect("equal")

# fig2, ax2 = plt.subplots()
# ax2.contourf(plotObj.gridX, plotObj.gridY, plotObj.gridVx)
# ax2.quiver(plotObj.gridX[::24,::24], plotObj.gridY[::24,::24],plotObj.gridVx[::24,::24],plotObj.gridVy[::24,::24])

# ax2.plot(discLoc,3.959,"kx")
# ax2.set_aspect("equal")

# fig3Obj = combinedFrameObj
# fig3, ax3 = plt.subplots()
# ax3.contourf(fig3Obj.gridX, fig3Obj.gridY, fig3Obj.gridWz, levels=np.linspace(-75,75,150),extend="both")
# ax3.quiver(fig3Obj.gridX[::24,::24], fig3Obj.gridY[::24,::24],fig3Obj.gridVx[::24,::24],fig3Obj.gridVy[::24,::24])
# ax3.set_aspect("equal")

# fig4, ax4 = plt.subplots()
# ax4.contourf(fig3Obj.gridX, fig3Obj.gridY, fig3Obj.computeWz(), levels=np.linspace(-75,75,150),extend="both")
# ax4.set_aspect("equal")

# fig5, ax5 = plt.subplots()
# ax5.contourf(fig3Obj.gridX, fig3Obj.gridY, curl2D(fig3Obj.gridX*1e-3,fig3Obj.gridY*1e-3,fig3Obj.gridVx,fig3Obj.gridVy), levels=np.linspace(-75,75,150),extend="both")
# ax5.set_aspect("equal")

fig4Obj = flipCombinedFrameObj
# fig6, ax6 = plt.subplots()
# ax6.contourf(fig4Obj.gridX, fig4Obj.gridY, fig4Obj.gridWz, levels=np.linspace(-75,75,150),extend="both")
# ax6.quiver(fig4Obj.gridX[::24,::24], fig4Obj.gridY[::24,::24],fig4Obj.gridVx[::24,::24],fig4Obj.gridVy[::24,::24])
# ax6.set_aspect("equal")

# %%

# fieldObj.getIndVel(0.05,0.05)
# indVxRaw1 = fieldObj.indVxRaw
# indVxIntp1 = fieldObj.indVxIntp

# fieldObj.getIndVel(0.05,0.1)
# indVxRaw2 = fieldObj.indVxRaw
# indVxIntp2 = fieldObj.indVxIntp

# fieldObj.getIndVel(0.05,0.15)
# indVxRaw3 = fieldObj.indVxRaw
# indVxIntp3 = fieldObj.indVxIntp

# fieldObj.getIndVel(0.05,0.2)
# indVxRaw4 = fieldObj.indVxRaw
# indVxIntp4 = fieldObj.indVxIntp

# fieldObj.getIndVel(0.05,0.25)
# indVxRaw5 = fieldObj.indVxRaw
# indVxIntp5 = fieldObj.indVxIntp

figInd, axInd = plt.subplots()
# axInd.plot(fieldObj.lineRaw[:,1],fieldObj.indVxRaw)
# axInd.plot(fieldObj.lineIntp[:,1],fieldObj.indVxIntp)
# axInd.plot(fieldObj.lineIntp[:,1],indVxIntp1)
# axInd.plot(fieldObj.lineIntp[:,1],indVxIntp2)
# axInd.plot(fieldObj.lineIntp[:,1],indVxIntp3)
# axInd.plot(fieldObj.lineIntp[:,1],indVxIntp4)
# axInd.plot(fieldObj.lineIntp[:,1],indVxIntp5)

#%%

# fieldObj.getIndVel(0.05,0.3)
# indVxRaw6 = fieldObj.indVxRaw
# indVxIntp6 = fieldObj.indVxIntp

# fieldObj.getIndVel(0.05,0.35)
# indVxRaw7 = fieldObj.indVxRaw
# indVxIntp7 = fieldObj.indVxIntp

# fieldObj.getIndVel(0.05,0.4)
# indVxRaw8 = fieldObj.indVxRaw
# indVxIntp8 = fieldObj.indVxIntp

# fieldObj.getIndVel(0.05,0.45)
# indVxRaw9 = fieldObj.indVxRaw
# indVxIntp9 = fieldObj.indVxIntp

# fieldObj.getIndVel(0.05,0.5)
# indVxRaw10 = fieldObj.indVxRaw
# indVxIntp10 = fieldObj.indVxIntp

# axInd.plot(fieldObj.lineIntp[:,1],indVxIntp6)
# axInd.plot(fieldObj.lineIntp[:,1],indVxIntp7)
# axInd.plot(fieldObj.lineIntp[:,1],indVxIntp8)
# axInd.plot(fieldObj.lineIntp[:,1],indVxIntp9)
# axInd.plot(fieldObj.lineIntp[:,1],indVxIntp10)

gridRndm = np.zeros(fig4Obj.gridX.shape)
gridRndm2 = np.zeros(fig4Obj.gridX.shape)
discMask = fieldObj.getDiscMask(fieldObj.field)
gridRndm[discMask[0]] = 1
gridRndm2[discMask[1]] = 2

fig7, ax7 = plt.subplots()
ax7.contourf(fig4Obj.gridX, fig4Obj.gridY,gridRndm)
ax7.contourf(fig4Obj.gridX, fig4Obj.gridY,gridRndm2, alpha = 0.5)

#%%

# figInd2, axInd2 = plt.subplots()
# # axInd.plot(fieldObj.lineRaw[:,1],fieldObj.indVxRaw)
# # axInd.plot(fieldObj.lineIntp[:,1],fieldObj.indVxIntp)
# axInd2.plot(fieldObj.lineIntp[:,1],indVxIntp1)
# axInd2.plot(fieldObj.lineIntp[:,1],indVxIntp2)
# axInd2.plot(fieldObj.lineIntp[:,1],indVxIntp3)
# axInd2.plot(fieldObj.lineIntp[:,1],indVxIntp4)
# axInd2.plot(fieldObj.lineIntp[:,1],indVxIntp5)
# axInd2.plot(fieldObj.lineIntp[:,1],indVxIntp6)
# axInd2.plot(fieldObj.lineIntp[:,1],indVxIntp7)
# axInd2.plot(fieldObj.lineIntp[:,1],indVxIntp8)
# axInd2.plot(fieldObj.lineIntp[:,1],indVxIntp9)
# axInd2.plot(fieldObj.lineIntp[:,1],indVxIntp10)

# rbfMask = fieldObj.getRBFMask(fieldObj.field)
rbfMask2 = fieldObj.getRBFMask2(fieldObj.field,ustrOffset=0.05,dstrOffset=0.25)

ustrMinY, ustrMaxY, ustrMinX, ustrMaxX = bool2DMinMax(rbfMask2[0])
dstrMinY, dstrMaxY, dstrMinX, dstrMaxX = bool2DMinMax(rbfMask2[1])

gridX = fieldObj.field.gridX
gridY = fieldObj.field.gridY

fig10,ax10 = plt.subplots()
ax10.contourf(gridX, gridY, gridRndm2)
ax10.contourf(gridX[ustrMinY:ustrMaxY,ustrMinX:ustrMaxX],gridY[ustrMinY:ustrMaxY,ustrMinX:ustrMaxX],rbfMask2[0][ustrMinY:ustrMaxY,ustrMinX:ustrMaxX])
ax10.contourf(gridX[dstrMinY:dstrMaxY,dstrMinX:dstrMaxX],gridY[dstrMinY:dstrMaxY,dstrMinX:dstrMaxX],rbfMask2[1][dstrMinY:dstrMaxY,dstrMinX:dstrMaxX])
ax10.set_aspect("equal")

# %%
# fieldObj.field.genRBFInterpolator(rbfMask, skipDims=4)
fieldObj.field.genRBFInterpolator2(*rbfMask2[:2], skipDims=3, rbfKernel="linear")

# %%
rawGridVx = fieldObj.field.gridVx
inpaintGridVx = np.array(rawGridVx)
# inpaintGridVx[discMask[0]] = fieldObj.field.gridVxIntpr(np.column_stack((gridX[discMask[0]],gridY[discMask[0]])))
inpaintGridVx[rbfMask2[3]] = fieldObj.field.gridVxRBFIntpr(np.column_stack((gridX[rbfMask2[3]],gridY[rbfMask2[3]])))

fig8, ax8 = plt.subplots()
ax8.contourf(gridX, gridY, rawGridVx)
ax8.axvline(flipDiscLoc, color="k", linestyle="--")
ax8.set_aspect("equal")

fig9, ax9 = plt.subplots()
ax9.contourf(gridX, gridY, inpaintGridVx)
ax9.axvline(flipDiscLoc, color="k", linestyle="--")
ax9.set_aspect("equal")
# %%
