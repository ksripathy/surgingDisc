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

from matplotlib import pyplot as plt
import numpy as np
from numpy.ma import masked_array
from mpl_toolkits.axes_grid1 import make_axes_locatable

from src.plotInd import plotInd
from src.ioUtilsNoStitch import caseData
from src.miscTools import axesBboxSize
from src.plotFields import plotFields

p45Case0 = plotInd(plotDir + "/noStitchRes/results_p45Case0.json")
p45Case1 = plotInd(plotDir + "/noStitchRes/results_p45Case1.json")
p45Case2 = plotInd(plotDir + "/noStitchRes/results_p45Case2.json")
p45Case3 = plotInd(plotDir + "/noStitchRes/results_p45Case3.json")
p45Case4 = plotInd(plotDir + "/noStitchRes/results_p45Case4.json")
p45Case5 = plotInd(plotDir + "/noStitchRes/results_p45Case5.json")
p45Case6 = plotInd(plotDir + "/noStitchRes/results_p45Case6.json")
p45Case7 = plotInd(plotDir + "/noStitchRes/results_p45Case7.json")


p70Case0 = plotInd(plotDir + "/noStitchRes/results_p70Case0.json")
p70Case1 = plotInd(plotDir + "/noStitchRes/results_p70Case1.json")
p70Case2 = plotInd(plotDir + "/noStitchRes/results_p70Case2.json")
p70Case3 = plotInd(plotDir + "/noStitchRes/results_p70Case3.json")
p70Case4 = plotInd(plotDir + "/noStitchRes/results_p70Case4.json")
p70Case5 = plotInd(plotDir + "/noStitchRes/results_p70Case5.json")
p70Case6 = plotInd(plotDir + "/noStitchRes/results_p70Case6.json")
p70Case7 = plotInd(plotDir + "/noStitchRes/results_p70Case7.json")

'''p45Case0.surgeCycle6Span(plotDir)
p45Case1.surgeCycle6Span(plotDir)
p45Case2.surgeCycle6Span(plotDir)
p45Case3.surgeCycle6Span(plotDir)
p45Case4.surgeCycle6Span(plotDir)
p45Case5.surgeCycle6Span(plotDir)

p70Case0.surgeCycle6Span(plotDir)
p70Case1.surgeCycle6Span(plotDir)
p70Case2.surgeCycle6Span(plotDir)
p70Case3.surgeCycle6Span(plotDir)
p70Case4.surgeCycle6Span(plotDir)
p70Case5.surgeCycle6Span(plotDir)'''

'''p45Case6.discPlane4Phase(plotDir+"/torquePaper")
p45Case7.discPlane4Phase(plotDir+"/torquePaper")
p70Case6.discPlane4Phase(plotDir+"/torquePaper")
p70Case7.discPlane4Phase(plotDir+"/torquePaper")'''

#Plotting induction for different scenarios
spanLocIndex = 18

label0 = f"k = {p45Case0.data.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case0.data.reducedSurgeAmplitude}"
label1 = f"k = {p45Case1.data.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case1.data.reducedSurgeAmplitude}"
label2 = f"k = {p45Case2.data.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case2.data.reducedSurgeAmplitude}"
label3 = f"k = {p45Case3.data.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case3.data.reducedSurgeAmplitude}"
label4 = f"k = {p45Case4.data.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case4.data.reducedSurgeAmplitude}"
label5 = f"k = {p45Case5.data.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case5.data.reducedSurgeAmplitude}"
label6 = f"k," + r"$x_{amp}/D$ = " + f"{p45Case6.data.reducedSurgeFrequency}, " + r"$V_\infty$ = " + f"{p45Case6.data.freestreamVelocity}" + r" $ms^{-1}$"
label7 = f"k," + r"$x_{amp}/D$ = " + f"{p45Case7.data.reducedSurgeFrequency}, " + r"$V_\infty$ = " + f"{p45Case7.data.freestreamVelocity}" + r" $ms^{-1}$"

'''fig, ax = plt.subplots(ncols=2)
fig.set_size_inches(26/2.54, 10/2.54)
p45Line0 = p45Case0.phases(spanLocIndex, ax[0], True)
p45Line0[0][0].set_label(label0)
p45Line0[1].set_color("k")
p45Line0[1].set_linestyle("-")
p45Line0[1].set_label(label6)
p45Line1 = p45Case1.phases(spanLocIndex, ax[0], False)
p45Line1[0][0].set_label(label1)
p45Line2 = p45Case2.phases(spanLocIndex, ax[0], False)
p45Line2[0][0].set_label(label2)
p45Line3 = p45Case3.phases(spanLocIndex, ax[0], True)
p45Line3[0][0].set_label(label3)
p45Line3[0][0].set_linestyle("--")
p45Line3[1].set_color("k")
p45Line3[1].set_label(label7)
p45Line4 = p45Case4.phases(spanLocIndex, ax[0], False)
p45Line4[0][0].set_label(label4)
p45Line4[0][0].set_linestyle("--")
p45Line5 = p45Case5.phases(spanLocIndex, ax[0], False)
p45Line5[0][0].set_label(label5)
p45Line5[0][0].set_linestyle("--")

ax[0].axvline(p45Case0.data.minLocTime, color = "r", linestyle = "dotted")
ax[0].axvline(p45Case0.data.meanLocTime1, color = "b", linestyle = "dotted")
ax[0].axvline(p45Case0.data.maxLocTime, color = "g", linestyle = "dotted")
ax[0].axvline(p45Case0.data.meanLocTime2, color = "b", linestyle = "dotted")

ax[0].set_xlabel(r"t/T [-]")
ax[0].set_ylabel(r"a [-]")
ax[0].grid()

plt.figlegend(loc="upper center", ncols=4, bbox_to_anchor=(0.5,1.06))

p70Line0 = p70Case0.phases(spanLocIndex, ax[1], True)
p70Line0[0][0].set_label(label0)
p70Line0[1].set_color("k")
p70Line0[1].set_linestyle("-")
p70Line0[1].set_label(label6)
p70Line1 = p70Case1.phases(spanLocIndex, ax[1], False)
p70Line1[0][0].set_label(label1)
p70Line2 = p70Case2.phases(spanLocIndex, ax[1], False)
p70Line2[0][0].set_label(label2)
p70Line3 = p70Case3.phases(spanLocIndex, ax[1], True)
p70Line3[0][0].set_label(label3)
p70Line3[0][0].set_linestyle("--")
p70Line3[1].set_color("k")
p70Line3[1].set_label(label7)
p70Line4 = p70Case4.phases(spanLocIndex, ax[1], False)
p70Line4[0][0].set_label(label4)
p70Line4[0][0].set_linestyle("--")
p70Line5 = p70Case5.phases(spanLocIndex, ax[1], False)
p70Line5[0][0].set_label(label5)
p70Line5[0][0].set_linestyle("--")

ax[1].axvline(p70Case0.data.minLocTime, color = "r", linestyle = "dotted")
ax[1].axvline(p70Case0.data.meanLocTime1, color = "b", linestyle = "dotted")
ax[1].axvline(p70Case0.data.maxLocTime, color = "g", linestyle = "dotted")
ax[1].axvline(p70Case0.data.meanLocTime2, color = "b", linestyle = "dotted")

ax[1].set_xlabel(r"t/T [-]")
ax[1].set_ylabel(r"a [-]")
ax[1].grid()

fig.savefig(plotDir+f"/combinedCases/spanIndex{spanLocIndex}.png", dpi=300, bbox_inches="tight", pad_inches=0)'''

'''#Torque paper reduced frequency legend labels
label0k = f"k = {p45Case0.data.reducedSurgeFrequency}"
label1k = f"k = {p45Case1.data.reducedSurgeFrequency}"
label2k = f"k = {p45Case2.data.reducedSurgeFrequency}"
label3k = f"k = {p45Case3.data.reducedSurgeFrequency}"
label4k = f"k = {p45Case4.data.reducedSurgeFrequency}"
label5k = f"k = {p45Case5.data.reducedSurgeFrequency}"
label6k = f"k = " + f"{p45Case6.data.reducedSurgeFrequency}"
#label7 = f"k," + r"$x_{amp}/D$ = " + f"{p45Case7.data.reducedSurgeFrequency}, " + r"$V_\infty$ = " + f"{p45Case7.data.freestreamVelocity}" + r" $ms^{-1}$"

#Torque paper reduced amplitude legend labels
label0a = r"$x_{amp}/D$ = " + f"{p45Case0.data.reducedSurgeAmplitude}"
label1a = r"$x_{amp}/D$ = " + f"{p45Case1.data.reducedSurgeAmplitude}"
label2a = r"$x_{amp}/D$ = " + f"{p45Case2.data.reducedSurgeAmplitude}"
label3a = r"$x_{amp}/D$ = " + f"{p45Case3.data.reducedSurgeAmplitude}"
label4a = r"$x_{amp}/D$ = " + f"{p45Case4.data.reducedSurgeAmplitude}"
label5a = r"$x_{amp}/D$ = " + f"{p45Case5.data.reducedSurgeAmplitude}"
label6a = r"$x_{amp}/D$ = " + f"{p45Case6.data.reducedSurgeFrequency}"
#label7 = f"k," + r"$x_{amp}/D$ = " + f"{p45Case7.data.reducedSurgeFrequency}, " + r"$V_\infty$ = " + f"{p45Case7.data.freestreamVelocity}" + r" $ms^{-1}$"

#Torque paper reduced frequency plots
fig2, ax2 = plt.subplots()
p70Line0 = p70Case0.phases(spanLocIndex, ax2, True)
p70Line0[0][0].set_label(label0k)
p70Line0[1].set_color("k")
p70Line0[1].set_linestyle("-")
p70Line0[1].set_label(label6k)
p70Line2 = p70Case2.phases(spanLocIndex, ax2, False)
p70Line2[0][0].set_label(label2k)
p70Line4 = p70Case4.phases(spanLocIndex, ax2, True)
p70Line4[0][0].set_label(label4k)
p70Line4[0][0].set_linestyle("--")
p70Line4[1].set_color("k")
ax2.axvline(p70Case0.data.minLocTime, color = "r", linestyle = "dotted")
ax2.axvline(p70Case0.data.meanLocTime1, color = "b", linestyle = "dotted")
ax2.axvline(p70Case0.data.maxLocTime, color = "g", linestyle = "dotted")
ax2.axvline(p70Case0.data.meanLocTime2, color = "b", linestyle = "dotted")
ax2.grid()
ax2.legend()
fig2.savefig(plotDir+"/torquePaper/P70redFreq.png", dpi=300, bbox_inches="tight", pad_inches=0)

#Torque paper reduced amplitude plots
fig3, ax3 = plt.subplots()
p70Line1 = p70Case1.phases(spanLocIndex, ax3, True)
p70Line1[0][0].set_label(label1a)
p70Line1[1].set_color("k")
p70Line1[1].set_linestyle("-")
p70Line1[1].set_label(label6a)
p70Line2 = p70Case2.phases(spanLocIndex, ax3, False)
p70Line2[0][0].set_label(label2a)
p70Line3 = p70Case3.phases(spanLocIndex, ax3, True)
p70Line3[0][0].set_label(label3a)
p70Line3[0][0].set_linestyle("--")
p70Line3[1].set_color("k")
ax3.axvline(p70Case0.data.minLocTime, color = "r", linestyle = "dotted")
ax3.axvline(p70Case0.data.meanLocTime1, color = "b", linestyle = "dotted")
ax3.axvline(p70Case0.data.maxLocTime, color = "g", linestyle = "dotted")
ax3.axvline(p70Case0.data.meanLocTime2, color = "b", linestyle = "dotted")
ax3.grid()
ax3.legend()
fig3.savefig(plotDir+"/torquePaper/P70redAmp.png", dpi=300, bbox_inches="tight", pad_inches=0)

#Torque paper vmax plots
fig4,ax4 = plt.subplots()
p70Line0 = p70Case0.phases(spanLocIndex, ax4, True)
p70Line0[0][0].set_label(r"$V_{max}$ = 0.25")
p70Line0[1].set_color("k")
p70Line0[1].set_linestyle("-")
p70Line0[1].set_label(r"$V_{max}$ = 0")
p70Line1 = p70Case1.phases(spanLocIndex, ax4, False)
p70Line1[0][0].set_label(r"$V_{max}$ = 0.25")
p70Line2 = p70Case2.phases(spanLocIndex, ax4, False)
p70Line2[0][0].set_label(r"$V_{max}$ = 0.5")
p70Line3 = p70Case3.phases(spanLocIndex, ax4, True)
p70Line3[0][0].set_label(r"$V_{max}$ = 0.75")
p70Line3[0][0].set_linestyle("--")
p70Line3[1].set_color("k")
p70Line4 = p70Case4.phases(spanLocIndex, ax4, False)
p70Line4[0][0].set_label(r"$V_{max}$ = 0.75")
p70Line4[0][0].set_linestyle("--")
p70Line5 = p70Case5.phases(spanLocIndex, ax4, False)
p70Line5[0][0].set_label(r"$V_{max}$ = 1.0125")
p70Line5[0][0].set_linestyle("--")

ax4.axvline(p70Case0.data.minLocTime, color = "r", linestyle = "dotted")
ax4.axvline(p70Case0.data.meanLocTime1, color = "b", linestyle = "dotted")
ax4.axvline(p70Case0.data.maxLocTime, color = "g", linestyle = "dotted")
ax4.axvline(p70Case0.data.meanLocTime2, color = "b", linestyle = "dotted")
ax4.legend()
ax4.set_xlabel(r"t/T [-]")
ax4.set_ylabel(r"a [-]")
ax4.grid()
fig4.savefig(plotDir+"/torquePaper/P70vmax.png", dpi=300, bbox_inches="tight", pad_inches=0)'''

'''#Flowfield plots
figData = np.empty(12, dtype=np.ndarray)
figsVel = np.empty(12, dtype=np.ndarray)
axsVel = np.empty(12, dtype=np.ndarray)
figsVort = np.empty(12, dtype=np.ndarray)
axsVort = np.empty(12, dtype=np.ndarray)

testCase = caseData(dataDir + "/tecPlotNoStitch/p70Case6*")

cdict = {'red': [[0.0, 0.0, 0.0],
                 [0.5, 1.0, 1.0],
                 [1.0, 0.0, 0.0]],
         'green': [[0.0, 0.0, 0.0],
                   [0.5, 1.0, 1.0],
                   [1.0, 0.0, 0.0]],
         'blue': [[0.0, 0.0, 0.0],
                  [0.5, 1.0, 1.0],
                  [1.0, 0.0, 0.0]],}

#newcmp = LinearSegmentedColormap('divGreys', segmentdata=cdict, N=125)

phases = np.arange(18, 360, 36)

for i in range(testCase.phaseQty):
    
    figData[i] = getattr(testCase,f"phase{i}").combinedFrames
    maskedX = masked_array(figData[i].gridPosX, mask=~figData[i].gridIsValid)
    maskedY = masked_array(figData[i].gridPosY, mask=~figData[i].gridIsValid)
    maskedVelX = masked_array(figData[i].gridVelX, mask=~figData[i].gridIsValid)
    maskedVelY = masked_array(figData[i].gridVelY, mask=~figData[i].gridIsValid)
    maskedVel = masked_array(figData[i].gridVel, mask=~figData[i].gridIsValid)
    maskedVort = masked_array(figData[i].gridVortZ, mask=~figData[i].gridIsValid)

    figsVel[i], axsVel[i] = plt.subplots()
    figObj = axsVel[i].contourf(figData[i].gridPosX/testCase.discDia, figData[i].gridPosY/testCase.discDia, maskedVel, levels=np.linspace(0.5,2.5,20), cmap="jet", extend="both")
    axsVel[i].quiver(maskedX[::24, ::24]/testCase.discDia, maskedY[::24, ::24]/testCase.discDia, maskedVelX[::24, ::24], maskedVelY[::24, ::24])
    divider2 = make_axes_locatable(axsVel[i])
    caxObj = divider2.append_axes("right", size="2%", pad=0.05)
    cbar = figsVel[i].colorbar(figObj, cax=caxObj)
    cbar.set_label(r"Velocity $\|V\|\ [ms^{-1}]$")
    axsVel[i].set_xlim(left=-350/testCase.discDia, right=300/testCase.discDia)
    axsVel[i].set_ylim(bottom =-170/testCase.discDia, top=180/testCase.discDia)
    axsVel[i].set_xlabel("x/D [-]")
    axsVel[i].set_ylabel("y/D [-]")
    axesBboxSize(6.5 + 0.02*6.5, 3.5, axsVel[i])
    
    figsVort[i], axsVort[i] = plt.subplots()
    figObj = axsVort[i].contourf(figData[i].gridPosX/testCase.discDia, figData[i].gridPosY/testCase.discDia, maskedVort, levels=np.linspace(-75,75,50), cmap="jet", extend="both")
    axsVort[i].quiver(maskedX[::24, ::24]/testCase.discDia, maskedY[::24, ::24]/testCase.discDia, maskedVelX[::24, ::24], maskedVelY[::24, ::24])
    divider2 = make_axes_locatable(axsVort[i])
    caxObj = divider2.append_axes("right", size="2%", pad=0.05)
    cbar = figsVort[i].colorbar(figObj, cax=caxObj)
    cbar.set_label(r"Vorticity $\omega_z\ [s^{-1}]$")
    axsVort[i].set_xlim(left=-350/testCase.discDia, right=300/testCase.discDia)
    axsVort[i].set_ylim(bottom =-170/testCase.discDia, top=180/testCase.discDia)
    axsVort[i].set_xlabel("x/D [-]")
    axsVort[i].set_ylabel("y/D [-]")
    axesBboxSize(6.5 + 0.02*6.5, 3.5, axsVort[i])

plt.show()'''

fieldPlots = plotFields(dataDir + "/tecPlotNoStitch/p70Case5*", axBounds = [-350, 300, 180, -170])

"Static case field plots"

'''staticVelFig, staticVelAx = plt.subplots(1, 2, figsize=(2 * 0.01 * 1.02 * (300 + 350), 1 * 0.01 * (180 + 170)), sharex=True, sharey=True, layout="constrained")
staticVortFig, staticVortAx = plt.subplots(1, 2, figsize=(2 * 0.01 * 1.02 * (300 + 350), 1 * 0.01 * (180 + 170)), sharex=True, sharey=True, layout="constrained")
case6FieldObj = plotFields(dataDir + "/tecPlotNoStitch/p70Case6*", axBounds = [-350, 300, 180, -170], figObjExt=[staticVelFig, staticVortFig])
case7FieldObj = plotFields(dataDir + "/tecPlotNoStitch/p70Case7*", axBounds = [-350, 300, 180, -170], figObjExt=[staticVelFig, staticVortFig])

cbarVel = staticVelFig.colorbar(case7FieldObj.velPlot, ax=staticVelAx[1], pad=0)
cbarVel.set_label(r"$\|V\|/V_\infty\ [-]$")

cbarVort = staticVortFig.colorbar(case7FieldObj.vortPlot, ax=staticVortAx[1], pad=0)
cbarVort.set_label(r"$\omega_z\ \times D/V_\infty [-]$")

staticVelAx[0].set_xlabel("x/D [-]")
staticVelAx[1].set_xlabel("x/D [-]")
staticVelAx[0].set_ylabel("y/D [-]")

staticVortAx[0].set_xlabel("x/D [-]")
staticVortAx[1].set_xlabel("x/D [-]")
staticVortAx[0].set_ylabel("y/D [-]")

staticVelAx[0].set_title(r"$V_\infty = 3ms^{-1}$")
staticVelAx[1].set_title(r"$V_\infty = 2ms^{-1}$")

staticVortAx[0].set_title(r"$V_\infty = 3ms^{-1}$")
staticVortAx[1].set_title(r"$V_\infty = 2ms^{-1}$")

staticVelFig.savefig(plotDir+"/fields/p70StaticVel.png", dpi=300)
staticVortFig.savefig(plotDir+"/fields/p70StaticVort.png", dpi=300)'''
