import os
import sys
#%matplotlib widget

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

#Add src folder to python path
#sys.path.append(srcDir)

import numpy as np
import matplotlib.pyplot as plt

from src.ioUtilsNoStitch import caseData

'''p45Case0 = caseData(dataDir + "/p45Case0*")
p45Case0.save(plotDir)
p45Case1 = caseData(dataDir + "/p45Case1*")
p45Case1.save(plotDir)
p45Case2 = caseData(dataDir + "/p45Case2*")
p45Case2.save(plotDir)
p45Case6 = caseData(dataDir + "/p45Case6*")
p45Case6.save(plotDir)

p45Case3 = caseData(dataDir + "/p45Case3*")
p45Case3.save(plotDir)
p45Case4 = caseData(dataDir + "/p45Case4*")
p45Case4.save(plotDir)
p45Case5 = caseData(dataDir + "/p45Case5*")
p45Case5.save(plotDir)

p70Case0 = caseData(dataDir + "/p70Case0*")
p70Case0.save(plotDir)
p70Case1 = caseData(dataDir + "/p70Case1*")
p70Case1.save(plotDir)
p70Case2 = caseData(dataDir + "/p70Case2*")
p70Case2.save(plotDir)
p70Case6 = caseData(dataDir + "/p70Case6*")
p70Case6.save(plotDir)

p70Case3 = caseData(dataDir + "/p70Case3*")
p70Case3.save(plotDir)'''
'''p70Case4 = caseData(dataDir + "/p70Case4*")
p70Case4.save(plotDir)'''
'''p70Case5 = caseData(dataDir + "/p70Case5*")
p70Case5.save(plotDir)

p45Case7 = caseData(dataDir + "/p45Case7*")
p45Case7.save(plotDir)

p70Case7 = caseData(dataDir + "/p70Case7*")
p70Case7.save(plotDir)'''

'''p45Case0 = caseData(dataDir + "/tecPlotNoStitch/p45Case0*")
p45Case0.saveAsMat(dataDir + "/matlab")
p45Case1 = caseData(dataDir + "/tecPlotNoStitch/p45Case1*")
p45Case1.saveAsMat(dataDir + "/matlab")
p45Case2 = caseData(dataDir + "/tecPlotNoStitch/p45Case2*")
p45Case2.saveAsMat(dataDir + "/matlab")
p45Case3 = caseData(dataDir + "/tecPlotNoStitch/p45Case3*")
p45Case3.saveAsMat(dataDir + "/matlab")
p45Case4 = caseData(dataDir + "/tecPlotNoStitch/p45Case4*")
p45Case4.saveAsMat(dataDir + "/matlab")
p45Case5 = caseData(dataDir + "/tecPlotNoStitch/p45Case5*")
p45Case5.saveAsMat(dataDir + "/matlab")
p45Case6 = caseData(dataDir + "/tecPlotNoStitch/p45Case6*")
p45Case6.saveAsMat(dataDir + "/matlab")
p45Case7 = caseData(dataDir + "/tecPlotNoStitch/p45Case7*")
p45Case7.saveAsMat(dataDir + "/matlab")'''

'''p70Case0 = caseData(dataDir + "/tecPlotNoStitch/p70Case0*")
p70Case0.saveAsMat(dataDir + "/matlab")
p70Case1 = caseData(dataDir + "/tecPlotNoStitch/p70Case1*")
p70Case1.saveAsMat(dataDir + "/matlab")
p70Case2 = caseData(dataDir + "/tecPlotNoStitch/p70Case2*")
p70Case2.saveAsMat(dataDir + "/matlab")
p70Case3 = caseData(dataDir + "/tecPlotNoStitch/p70Case3*")
p70Case3.saveAsMat(dataDir + "/matlab")'''
p70Case4 = caseData(dataDir + "/tecPlotNoStitch/p70Case4*")
p70Case4.saveAsMat(dataDir + "/matlab")
p70Case5 = caseData(dataDir + "/tecPlotNoStitch/p70Case5*")
p70Case5.saveAsMat(dataDir + "/matlab")
'''p70Case6 = caseData(dataDir + "/tecPlotNoStitch/p70Case6*")
p70Case6.saveAsMat(dataDir + "/matlab")
p70Case7 = caseData(dataDir + "/tecPlotNoStitch/p70Case7*")
p70Case7.saveAsMat(dataDir + "/matlab")'''

