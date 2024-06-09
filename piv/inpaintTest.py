import os
import sys
import numpy as np
from PIL import Image
import scipy
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

from piv.src.pivFields import planarPIVField

p70Case7 = planarPIVField(dataDir+"/tecPlotNoStitch/p70Case7Mean0001_0.dat")
p70Case7.combine2Frames()

testNan = np.array(p70Case7.combinedFrames.gridVel)
testNan[~p70Case7.combinedFrames.gridIsValid] = np.nan

ipn_kernel = np.array([[1,1,1],[1,0,1],[1,1,1]]) # kernel for inpaint_nans

def inpaint_nans(im):
    nans = np.isnan(im)
    while np.sum(nans)>0:
        im[nans] = 0
        vNeighbors = scipy.signal.convolve2d((nans==False),ipn_kernel,mode='same',boundary='symm')
        im2 = scipy.signal.convolve2d(im,ipn_kernel,mode='same',boundary='symm')
        im2[vNeighbors>0] = im2[vNeighbors>0]/vNeighbors[vNeighbors>0]
        im2[vNeighbors==0] = np.nan
        im2[(nans==False)] = im[(nans==False)]
        im = im2
        nans = np.isnan(im)
    return im

testNaNFilled = inpaint_nans(testNan)

x = p70Case7.frameAxialLocs
y = p70Case7.frameSpanLocs
xx,yy = np.meshgrid(x,y)
x1 = xx[p70Case7.combinedFrames.gridIsValid]
y1 = yy[p70Case7.combinedFrames.gridIsValid]
newarr = p70Case7.combinedFrames.gridVel[p70Case7.combinedFrames.gridIsValid]
GD1 = scipy.interpolate.griddata((x1,y1), newarr.ravel(), (xx, yy), method="cubic")

gridVelFilled = np.array(p70Case7.combinedFrames.gridVel)
gridVelFilled[~p70Case7.combinedFrames.gridIsValid] = scipy.interpolate.griddata((x1,y1), newarr.ravel(), (xx[~p70Case7.combinedFrames.gridIsValid], yy[~p70Case7.combinedFrames.gridIsValid]), method="cubic")

fig1, ax1 = plt.subplots()
ax1.contourf(p70Case7.combinedFrames.gridPosX,p70Case7.combinedFrames.gridPosY,p70Case7.combinedFrames.gridVel)

fig2, ax2 = plt.subplots()
ax2.contourf(p70Case7.combinedFrames.gridPosX,p70Case7.combinedFrames.gridPosY,testNaNFilled)

fig4, ax4 = plt.subplots()
ax4.contourf(p70Case7.combinedFrames.gridPosX,p70Case7.combinedFrames.gridPosY,gridVelFilled)

p70Case7.fillMask(-1)

fig5, ax5 = plt.subplots()
ax5.contourf(p70Case7.combinedFrames.gridPosX,p70Case7.combinedFrames.gridPosY,p70Case7.combinedFrames.gridVelFilled, cmap="seismic")

'''fig3, ax3 = plt.subplots()
ax3.imshow(~p70Case6.combinedFrames.gridIsValid.astype(np.uint8))

img = Image.fromarray(~p70Case6.combinedFrames.gridIsValid,"1")
img.show()'''
