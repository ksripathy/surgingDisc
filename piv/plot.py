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
from src.plotInd import plotInd

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

p45Case6.discPlane4Phase(plotDir)
p45Case7.discPlane4Phase(plotDir)
p70Case6.discPlane4Phase(plotDir)
p70Case7.discPlane4Phase(plotDir)

'''#Plotting induction for different scenarios
spanLocIndex = 18

label0 = f"k = {p45Case0.data.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case0.data.reducedSurgeAmplitude}"
label1 = f"k = {p45Case1.data.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case1.data.reducedSurgeAmplitude}"
label2 = f"k = {p45Case2.data.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case2.data.reducedSurgeAmplitude}"
label3 = f"k = {p45Case3.data.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case3.data.reducedSurgeAmplitude}"
label4 = f"k = {p45Case4.data.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case4.data.reducedSurgeAmplitude}"
label5 = f"k = {p45Case5.data.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{p45Case5.data.reducedSurgeAmplitude}"
label6 = f"k," + r"$x_{amp}/D$ = " + f"{p45Case6.data.reducedSurgeFrequency}, " + r"$V_\infty$ = " + f"{p45Case6.data.freestreamVelocity}" + r" $ms^{-1}$"
label7 = f"k," + r"$x_{amp}/D$ = " + f"{p45Case7.data.reducedSurgeFrequency}, " + r"$V_\infty$ = " + f"{p45Case7.data.freestreamVelocity}" + r" $ms^{-1}$"

fig, ax = plt.subplots(ncols=2)
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



