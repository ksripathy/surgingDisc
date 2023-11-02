import os
import sys
#%matplotlib widget

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

#Add src folder to python path
sys.path.append(srcDir)

import numpy as np
from numpy.ma import masked_array
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import UnivariateSpline

from src.pivFields import planarPIVField
from src.miscTools import axesBboxSize
from src.miscTools import attributeArray
from src.miscTools import nearestValueIndex
from src.ioUtilsNoStitch import caseData
from src.ioUtilsNoStitch import loadData

'''#p45Case6 = caseData(dataDir + "/tecPlotNoStitch/p45Case0*")
#p45Case6.save(plotDir + "/noStitchRes")
p45Case6 = loadData(plotDir + "/noStitchRes/results_p45Case0.json")
fig5Data = p45Case6.cycleAxialIndV2
fig6Data = p45Case6.cycleAxialInd
phase = 2

fig5, ax5 = plt.subplots()
ax5.plot(fig5Data["spanLocs"], fig5Data["intp0"][phase], label="intp0", linestyle="--", color="k")
#ax5.plot(fig5Data["spanLocs"], fig5Data["intp5"][phase], label="intp5")
#ax5.plot(fig5Data["spanLocs"], fig5Data["intp25"][phase], label="intp25")
#ax5.plot(fig5Data["spanLocs"], fig5Data["intp50"][phase], label="intp50")
ax5.plot(fig5Data["spanLocs"], fig5Data["intp100"][phase], label="intp100")
ax5.plot(fig5Data["spanLocs"], fig5Data["intp200"][phase], label="intp200")
plt.legend()

fig6, ax6 = plt.subplots()
ax6.plot(fig6Data["spanLocs"], fig6Data["intp0"][phase], label="intp0", linestyle="--", color="k")
#ax6.plot(fig6Data["spanLocs"], fig6Data["intp5"][phase], label="intp5")
#ax6.plot(fig6Data["spanLocs"], fig6Data["intp25"][phase], label="intp25")
#ax6.plot(fig6Data["spanLocs"], fig6Data["intp50"][phase], label="intp50")
ax6.plot(fig6Data["spanLocs"], fig6Data["intp100"][phase], label="intp100")
ax6.plot(fig6Data["spanLocs"], fig6Data["intp200"][phase], label="intp200")
plt.legend()'''


testCase = planarPIVField(dataDir + "/tecPlotNoStitch/p45Case2Mean0001_6.dat")
#testCase2 = planarPIVField(dataDir + "/p45Case4Mean0001_2.dat")
testCase.combine2Frames(overlapSmoothing=True)

spanLoc = 5.962
spanLocs = (np.linspace(spanLoc-100, spanLoc+100, 21) - spanLoc)/200
spanLocIndex = nearestValueIndex(testCase.frame1.gridPosY[:,0], spanLoc)

maskBoundIndex = testCase.maskBoundaryIndex(1,spanLocIndex, calbEdgeArtifactTrim=5)
#axialLoc = testCase.frame1.gridPosX[spanLocIndex, maskBoundIndex]
axialLoc = 0

'''axialIndDistr3v2 = testCase.discAxialIndDistrV2(1, spanLoc, axialLoc, 3)
axialIndDistr5v2 = testCase.discAxialIndDistrV2(1, spanLoc, axialLoc, 5)
axialIndDistr10v2 = testCase.discAxialIndDistrV2(1, spanLoc, axialLoc, 10)
axialIndDistr20v2 = testCase.discAxialIndDistrV2(1, spanLoc, axialLoc, 20)
axialIndDistr40v2 = testCase.discAxialIndDistrV2(1, spanLoc, axialLoc, 40)
axialIndDistr0v2 = testCase.discAxialIndDistrV2(1, spanLoc, axialLoc, 0)

fig3, ax3 = plt.subplots()
ax3.plot(spanLocs, axialIndDistr3v2, label="intp3")
ax3.plot(spanLocs, axialIndDistr5v2, label="intp5")
ax3.plot(spanLocs, axialIndDistr10v2, label="intp10")
ax3.plot(spanLocs, axialIndDistr20v2, label="intp20")
ax3.plot(spanLocs, axialIndDistr40v2, label="intp40")
ax3.plot(spanLocs, axialIndDistr0v2, label="intp0", linestyle="--", color="k")
ax3.legend()

axialIndDistr3v1 = testCase.discAxialIndDistr(1, spanLoc, axialLoc, 3)
axialIndDistr5v1 = testCase.discAxialIndDistr(1, spanLoc, axialLoc, 5)
axialIndDistr10v1 = testCase.discAxialIndDistr(1, spanLoc, axialLoc, 10)
axialIndDistr20v1 = testCase.discAxialIndDistr(1, spanLoc, axialLoc, 20)
axialIndDistr40v1 = testCase.discAxialIndDistr(1, spanLoc, axialLoc, 40)
axialIndDistr0v1 = testCase.discAxialIndDistr(1, spanLoc, axialLoc, 0)

fig4, ax4 = plt.subplots()
ax4.plot(spanLocs, axialIndDistr3v1, label="intp3")
ax4.plot(spanLocs, axialIndDistr5v1, label="intp5")
ax4.plot(spanLocs, axialIndDistr10v1, label="intp10")
ax4.plot(spanLocs, axialIndDistr20v1, label="intp20")
ax4.plot(spanLocs, axialIndDistr40v1, label="intp40")
ax4.plot(spanLocs, axialIndDistr0v1, label="intp0", linestyle="--", color="k")
ax4.legend()'''

axialLocQty =200

fig7X = np.linspace(axialLoc, axialLoc+axialLocQty-1, axialLocQty)
axialIndAxialDistr = []
axialIndPolyFit2 = []
axialIndPolyFit3 = []
axialIndPolyFit4 = []
axialIndPolyFit5 = []
axialLocMaskEnd = testCase.frame1.gridPosX[0,maskBoundIndex ]
for i in range(axialLocQty):
    
    temp = testCase.axialIndSpanDistr(1, spanLoc, axialLocMaskEnd + i)
    
    for j in range(len(temp)):
        
        if i == 0:
            
            axialIndAxialDistr.append([temp[j]])
            
        else:
            
            axialIndAxialDistr[j].append(temp[j])
            
        if i == axialLocQty-1:
            
            #axialIndPolyFit2.append(np.polyfit(fig7X, axialIndAxialDistr[j], 2))
            axialIndPolyFit3.append(np.polyfit(fig7X, axialIndAxialDistr[j], 3))
            axialIndPolyFit4.append(np.polyfit(fig7X, axialIndAxialDistr[j], 4))
            #axialIndPolyFit5.append(np.polyfit(fig7X, axialIndAxialDistr[j], 5))
            

'''fig7SpanIndex = 0
    
fig7, ax7 = plt.subplots()
ax7.plot(fig7X, axialIndAxialDistr[fig7SpanIndex])
ax7.plot(fig7X, np.poly1d(axialIndPolyFit2[fig7SpanIndex])(fig7X), label="k=2")
ax7.plot(fig7X, np.poly1d(axialIndPolyFit3[fig7SpanIndex])(fig7X), label="k=3")
#ax7.plot(fig7X, np.poly1d(axialIndPolyFit4[fig7SpanIndex])(fig7X), label="k=4")
#ax7.plot(fig7X, np.poly1d(axialIndPolyFit5[fig7SpanIndex])(fig7X), label="k=5")
#ax7.plot(fig7X, axialIndSplineFit[fig7SpanIndex](fig7X), label="spline")
ax7.legend()
for i in range(21):
    ax7.plot(fig7X, axialIndAxialDistr[i])'''
    
figs8 = []
for i in range(21):
    
    figs8.append(plt.figure())
    ax = figs8[i].add_axes(rect=(0.05,0.05,0.9,0.9))
    
    ax.plot(fig7X, axialIndAxialDistr[i])
    #ax.plot(fig7X, np.poly1d(axialIndPolyFit2[i])(fig7X), label="k=2")
    ax.plot(fig7X, np.poly1d(axialIndPolyFit3[i])(fig7X), label="k=3")
    ax.plot(fig7X, np.poly1d(axialIndPolyFit4[i])(fig7X), label="k=4")
    ax.legend()
    
    
    
    



maskCdict = {'red': [[0.0, 1.0, 1.0],
                     [1.0, 1.0, 1.0]],
         'green': [[0.0, 1.0, 1.0],
                   [1.0, 1.0, 1.0]],
         'blue': [[0.0, 1.0, 1.0],
                  [1.0, 1.0, 1.0]],}

maskCmp = LinearSegmentedColormap('maskWhite', segmentdata=maskCdict, N=125)

'''fig0,ax0 = plt.subplots()
#Masking areas with mask and other camera frame
fig0Data = masked_array(testCase.frame0.gridVortZ, mask=~testCase.frame0.gridIsValid)
fig0Plot = ax0.contourf(testCase.frame0.gridPosX, testCase.frame0.gridPosY, fig0Data, levels=np.linspace(-125,125,250), cmap="jet", extend="both")
#ax0.quiver(testCase.frame0.gridPosX[::16, ::16], testCase.frame0.gridPosY[::16, ::16], testCase.frame0.gridVelX[::16, ::16], testCase.frame0.gridVelY[::16, ::16])
divider0 = make_axes_locatable(ax0)
cax0 = divider0.append_axes("right", size="2%", pad=0.05) 
fig0.colorbar(fig0Plot, cax=cax0)
axesBboxSize(4.675 + 0.02*4.675, 2.237, ax0)
fig0.savefig("fig0.png", dpi=300, bbox_inches="tight")'''

'''fig1,ax1 = plt.subplots()
#Masking areas with mask and other camera frame
fig1Data = masked_array(testCase.frame1.gridVelX, mask=~testCase.frame1.gridIsValid)
fig1Plot = ax1.contourf(testCase.frame1.gridPosX, testCase.frame1.gridPosY, fig1Data, levels=100, cmap="jet", extend="both")
#ax1.quiver(testCase.frame1.gridPosX[::16, ::16], testCase.frame1.gridPosY[::16, ::16], testCase.frame1.gridVelX[::16, ::16], testCase.frame1.gridVelY[::16, ::16])
divider1 = make_axes_locatable(ax1)
cax1 = divider1.append_axes("right", size="2%", pad=0.05) 
fig1.colorbar(fig1Plot, cax=cax1)
axesBboxSize(4.675 + 0.02*4.675, 2.237, ax1)
fig1.savefig("fig1.png", dpi=300, bbox_inches="tight")'''

'''fig2,ax2 = plt.subplots()
#Masking areas with mask
fig2Data = masked_array(testCase.combinedFrames.gridVortZ, mask=~testCase.combinedFrames.gridIsValid)
fig2Plot = ax2.contourf(testCase.combinedFrames.gridPosX, testCase.combinedFrames.gridPosY, fig2Data, levels=np.linspace(-125,125,250), cmap="jet", extend="both")
#ax2.quiver(testCase.combinedFrames.gridPosX[::16, ::16], testCase.combinedFrames.gridPosY[::16, ::16], testCase.combinedFrames.gridVelX[::16, ::16], testCase.combinedFrames.gridVelY[::16, ::16])
divider2 = make_axes_locatable(ax2)
cax2 = divider2.append_axes("right", size="2%", pad=0.05) 
fig2.colorbar(fig2Plot, cax=cax2)
axesBboxSize(4.675 + 0.02*4.675, 2.237, ax2)
fig2.savefig("fig2.png", dpi=300, bbox_inches="tight")'''

plt.show()