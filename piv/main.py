import os
import sys
#%matplotlib widget

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots/noStitchRes")
dataDir = os.path.join(homeDir,"data")

#Add src folder to python path
sys.path.append(srcDir)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

from src.ioUtils import caseData
from src.ioUtils import loadData

#p45Case4 = caseData(dataDir + "/p45Case4*")

p45Case0 = loadData(plotDir + "/results_p45Case0.json")
p45Case1 = loadData(plotDir + "/results_p45Case1.json")
p45Case2 = loadData(plotDir + "/results_p45Case2.json")
p45Case3 = loadData(plotDir + "/results_p45Case3.json")
p45Case4 = loadData(plotDir + "/results_p45Case4.json")
p45Case5 = loadData(plotDir + "/results_p45Case5.json")
p45Case6 = loadData(plotDir + "/results_p45Case6.json")
p45Case7 = loadData(plotDir + "/results_p45Case7.json")

p70Case0 = loadData(plotDir + "/results_p70Case0.json")
p70Case1 = loadData(plotDir + "/results_p70Case1.json")
p70Case2 = loadData(plotDir + "/results_p70Case2.json")
p70Case3 = loadData(plotDir + "/results_p70Case3.json")
p70Case4 = loadData(plotDir + "/results_p70Case4.json")
p70Case5 = loadData(plotDir + "/results_p70Case5.json")
p70Case6 = loadData(plotDir + "/results_p70Case6.json")
p70Case7 = loadData(plotDir + "/results_p70Case7.json")

fig2, ax2 = plt.subplots(ncols=2)
fig2.set_size_inches(26/2.54, 10/2.54)
line1 = ax2[0].plot(p45Case0.ndimTimes, p45Case0.axialInductionDiscCentre, label=f"k = {p45Case0.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case0.reducedSurgeAmplitude}", marker = "x", color = "r")
line2 = ax2[0].plot(p45Case1.ndimTimes, p45Case1.axialInductionDiscCentre, label=f"k = {p45Case1.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case1.reducedSurgeAmplitude}", marker = "x", color = "g")
line3 = ax2[0].plot(p45Case2.ndimTimes, p45Case2.axialInductionDiscCentre, label=f"k = {p45Case2.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case2.reducedSurgeAmplitude}", marker = "x", color = "b")
line4 = ax2[0].plot(p45Case3.ndimTimes, p45Case3.axialInductionDiscCentre, label=f"k = {p45Case3.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case3.reducedSurgeAmplitude}", marker = "x", color = "c", linestyle = "--")
line5 = ax2[0].plot(p45Case4.ndimTimes, p45Case4.axialInductionDiscCentre, label=f"k = {p45Case4.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case4.reducedSurgeAmplitude}", marker = "x", color = "m", linestyle = "--")
line6 = ax2[0].plot(p45Case5.ndimTimes, p45Case5.axialInductionDiscCentre, label=f"k = {p45Case5.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case5.reducedSurgeAmplitude}", marker = "x", color = "y", linestyle = "--")
ax2[0].axhline(p45Case6.axialInductionDiscCentre, label=f"k," + r"$x_{amp}/D$ = " + f"{p45Case6.reducedSurgeFrequency}, " + r"$V_\infty$ = " + f"{p45Case6.freestreamVelocity}" + r" $ms^{-1}$", color = "k")
ax2[0].axhline(p45Case7.axialInductionDiscCentre, label=f"k," + r"$x_{amp}/D$ = " + f"{p45Case7.reducedSurgeFrequency}, " + r"$V_\infty$ = " + f"{p45Case7.freestreamVelocity}" + r" $ms^{-1}$", color = "k", linestyle = "--")
ax2[0].axvline(p45Case0.minLocTime, color = "r", linestyle = "dotted")
ax2[0].axvline(p45Case0.meanLocTime1, color = "b", linestyle = "dotted")
ax2[0].axvline(p45Case0.maxLocTime, color = "g", linestyle = "dotted")
ax2[0].axvline(p45Case0.meanLocTime2, color = "b", linestyle = "dotted")
#ax2.legend(loc="upper center", bbox_to_anchor=(0.5,0))
ax2[0].grid()
#ax2[0].set_xlim(0,1)
ax2[0].set_xlabel(r"t/T [-]")
ax2[0].set_ylabel(r"a [-]")
plt.figlegend(loc="upper center", ncols=4, bbox_to_anchor=(0.5,1.05))
#fig2.savefig(plotDir + "/p45inductionRevised.png", dpi=300, bbox_inches="tight", pad_inches=0)

ax2[1].plot(p70Case0.ndimTimes, p70Case0.axialInductionDiscCentre, marker = "x", color = "r")
ax2[1].plot(p70Case1.ndimTimes, p70Case1.axialInductionDiscCentre, marker = "x", color = "g")
ax2[1].plot(p70Case2.ndimTimes, p70Case2.axialInductionDiscCentre, marker = "x", color = "b")
ax2[1].plot(p70Case3.ndimTimes, p70Case3.axialInductionDiscCentre, marker = "x", color = "c", linestyle = "--")
ax2[1].plot(p70Case4.ndimTimes, p70Case4.axialInductionDiscCentre, marker = "x", color = "m", linestyle = "--")
ax2[1].plot(p70Case5.ndimTimes, p70Case5.axialInductionDiscCentre, marker = "x", color = "y", linestyle = "--")
ax2[1].axhline(p70Case6.axialInductionDiscCentre, color = "k")
ax2[1].axhline(p70Case7.axialInductionDiscCentre, color = "k", linestyle = "--")
ax2[1].axvline(p70Case0.minLocTime, color = "r", linestyle = "dotted")
ax2[1].axvline(p70Case0.meanLocTime1, color = "b", linestyle = "dotted")
ax2[1].axvline(p70Case0.maxLocTime, color = "g", linestyle = "dotted")
ax2[1].axvline(p70Case0.meanLocTime2, color = "b", linestyle = "dotted")
#ax2[1].legend()
ax2[1].grid()
#ax2[1].set_xlim(0,1)
ax2[0].set_xlabel(r"t/T [-]")
ax2[1].set_ylabel(r"a [-]")
#fig3.savefig(plotDir + "/p70induction.png", dpi=300, bbox_inches="tight", pad_inches=0)
fig2.savefig(plotDir + "/p45p70induction.png", dpi=300, bbox_inches="tight", pad_inches=0)
testCase = caseData(dataDir + "/p45Case2*")
#axialLocs60, velXDistr60, velYDistr60 = testCase.phase0.velDistrAxial(60.0, 1500)
'''axialLocs72, velXDistr72, velYDistr72 = testCase.phase0.velDistrAxial(72.0, 1500)
#axialLocs84, velXDistr84, velYDistr84 = testCase.phase0.velDistrAxial(84.0, 1500)
fig1, ax1 = plt.subplots()
#ax1.plot(axialLocs60, velXDistr60)
ax1.plot(axialLocs72, velXDistr72)
#ax1.plot(axialLocs84, velXDistr84)
plt.show()'''

fig3, ax3 = plt.subplots()
ax3.plot(p70Case0.ndimTimes, p70Case0.axialInductionDiscCentre, marker = "x", color = "r", label=f"k = {p70Case0.reducedSurgeFrequency}")
ax3.plot(p70Case2.ndimTimes, p70Case2.axialInductionDiscCentre, marker = "x", color = "g", label=f"k = {p70Case2.reducedSurgeFrequency}")
ax3.plot(p70Case4.ndimTimes, p70Case4.axialInductionDiscCentre, marker = "x", color = "b", label=f"k = {p70Case4.reducedSurgeFrequency}", linestyle = "--")
ax3.axhline(p70Case6.axialInductionDiscCentre, color = "k", label=f"k =" + f"{p70Case6.reducedSurgeFrequency}")
ax3.axhline(p70Case7.axialInductionDiscCentre, color = "k", linestyle = "--", label=f"k =" + f"{p70Case7.reducedSurgeFrequency}")
ax3.legend()
ax3.grid()
ax3.set_xlabel(r"t/T [-]")
ax3.set_ylabel(r"a [-]")

#discLoc60, discVel60 = testCase.phase6.velDisc(60.0)
#discLoc72, discVel72 = testCase.phase6.velDisc(72.0, discLoc60)
#discLoc84, discVel84 = testCase.phase6.velDisc(84.0)

#ind60, ind72, ind84 = testCase.phase6.discAxialInd()

'''figData = np.empty(12, dtype=np.ndarray)
figs = np.empty(12, dtype=np.ndarray)
axs = np.empty(12, dtype=np.ndarray)

cdict = {'red': [[0.0, 0.0, 0.0],
                 [0.5, 1.0, 1.0],
                 [1.0, 0.0, 0.0]],
         'green': [[0.0, 0.0, 0.0],
                   [0.5, 1.0, 1.0],
                   [1.0, 0.0, 0.0]],
         'blue': [[0.0, 0.0, 0.0],
                  [0.5, 1.0, 1.0],
                  [1.0, 0.0, 0.0]],}

newcmp = LinearSegmentedColormap('divGreys', segmentdata=cdict, N=125)

for i in range(12):
    
    figData[i] = getattr(testCase,f"phase{i}")

    figs[i], axs[i] = plt.subplots()
    axs[i].contourf(figData[i].gridPosX, figData[i].gridPosY, figData[i].gridVortZ, levels=np.linspace(-125,125,4), cmap=newcmp, extend="both")
    #axs[i].set_title(f"phase{i}")
    axs[i].xaxis.set_visible(False)
    axs[i].yaxis.set_visible(False)
    figs[i].savefig(plotDir + f"/posterPhase{i}.png", dpi=300, bbox_inches="tight", pad_inches=0)

plt.show()'''


