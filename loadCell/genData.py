import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from src.miscTools import harmonicPassFilt
from src.miscTools import freqsSignal
from src.miscTools import cycleAvgStd

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

from src.ioUtilsRev import caseData

p45Case0 = caseData("45","02","20N",dataDir,1)
p45Case0.genAeroPlots()

# dftFreqs, surgeSPD = freqsSignal(p45Case0.inertialLoadsObj.dataExtdClean[:,0] - np.mean(p45Case0.inertialLoadsObj.dataExtdClean[:,0]),2e3)
# surgeFreq = dftFreqs[np.argmax(np.abs(surgeSPD))]



# cycleNo = 100
# plotAttrb1 = getattr(p45Case0.aeroLoadsObj,"cycleNdimIntpData")
# plotAttrb2 = getattr(p45Case0.aeroLoadsObj,"cycleNdimDirectFiltData")
# plotAttrb3 = getattr(p45Case0.aeroLoadsObj,"cycleNdimFiltIntpData")
# plotAttrb4 = getattr(p45Case0.aeroLoadsObj,"cycleNdimIntpFiltData")
# plotAttrb5 = getattr(p45Case0.aeroLoadsObj,"intpFiltCycleData")
# fig0, ax0 = plt.subplots()
# ax0.plot(plotAttrb1[cycleNo,2],plotAttrb1[cycleNo,4])
# ax0.plot(plotAttrb2[cycleNo,2],np.real(plotAttrb2[cycleNo,4].astype(complex)))
# ax0.plot(plotAttrb3[cycleNo,2],np.real(plotAttrb3[cycleNo,4].astype(complex)))
# ax0.plot(plotAttrb4[cycleNo,2],np.real(plotAttrb4[cycleNo,4].astype(complex)))
# ax0.plot(plotAttrb4[cycleNo,2],np.real(plotAttrb4[cycleNo,4].astype(complex)))
# ax0.plot(plotAttrb5[cycleNo,2],plotAttrb5[cycleNo,4])

# plotCase1P45 = caseData("45","18","2N",dataDir)
# plotCase2P45 = caseData("45","19","2N",dataDir)
# plotCase3P45 = caseData("45","20","2N",dataDir)
# plotCase4P45 = caseData("45","21","2N",dataDir)

# p45Cases = [plotCase1P45,plotCase2P45,plotCase3P45,plotCase4P45]

# plotCase1P70 = caseData("70","18","2N",dataDir)
# plotCase2P70 = caseData("70","19","2N",dataDir)
# plotCase3P70 = caseData("70","20","2N",dataDir)
# plotCase4P70 = caseData("70","21","2N",dataDir)

# p70Cases = [plotCase1P70,plotCase2P70,plotCase3P70,plotCase4P70]

# p45Fig, p45Ax = plt.subplots()
# p70Fig, p70Ax = plt.subplots()

# for p45Case, p70Case in zip(p45Cases,p70Cases):

#     p45Ax.plot(p45Case.aeroLoadsObj.ndimDirectFiltCycleAvgData[0],p45Case.aeroLoadsObj.ndimDirectFiltCycleAvgData[1],label=f"p{p45Case.discPor}Case{p45Case.caseNo}")
#     p45Ax.fill_between(p45Case.aeroLoadsObj.ndimDirectFiltCycleAvgData[0],p45Case.aeroLoadsObj.ndimDirectFiltCycleAvgData[1] - p45Case.aeroLoadsObj.ndimDirectFiltCycleStdData[1], p45Case.aeroLoadsObj.ndimDirectFiltCycleAvgData[1] + p45Case.aeroLoadsObj.ndimDirectFiltCycleStdData[1],alpha=0.5)
    
#     p70Ax.plot(p70Case.aeroLoadsObj.ndimDirectFiltCycleAvgData[0],p70Case.aeroLoadsObj.ndimDirectFiltCycleAvgData[1],label=f"p{p70Case.discPor}Case{p70Case.caseNo}")
#     p70Ax.fill_between(p70Case.aeroLoadsObj.ndimDirectFiltCycleAvgData[0],p70Case.aeroLoadsObj.ndimDirectFiltCycleAvgData[1] - p70Case.aeroLoadsObj.ndimDirectFiltCycleStdData[1], p70Case.aeroLoadsObj.ndimDirectFiltCycleAvgData[1] + p70Case.aeroLoadsObj.ndimDirectFiltCycleStdData[1],alpha=0.5)
    
# p45Ax.axhline(0.84,color="k",linestyle="--",label="p45Static")
# p45Ax.legend()
# p45Ax.set_xlabel(r"$\phi [\degree]$")
# p45Ax.set_ylabel(r"$C_T [-]$")

# p70Ax.axhline(0.46,color="k",linestyle="--",label="p70Static")
# p70Ax.legend()
# p70Ax.set_xlabel(r"$\phi [\degree]$")
# p70Ax.set_ylabel(r"$C_T [-]$")
