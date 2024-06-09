import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt

from src.recordedPhasesForcesRev import recordedData
from src.miscTools import dftPlot

class caseData:
    
    def __init__(self,discPor,caseNo,loadCell,dataDir,maxHarmonic=1):
        
        self.maxHarmonic = maxHarmonic
        self.caseNo = caseNo
        self.discPor = discPor
        
        caseParams = np.loadtxt(dataDir + "/surgeCaseParams.csv", delimiter=",")
        
        self.surgeAmp = caseParams[int(caseNo),3] #A
        self.surgeFreq = caseParams[int(caseNo),4] #f
        self.surgeAmpRed = caseParams[int(caseNo),0] #A*
        self.surgeFreqRed = caseParams[int(caseNo),1] #k
        self.surgeVelRed = caseParams[int(caseNo),2] #u*
        
        logFilePath = dataDir + f"/testMatrixTotal/curved/{loadCell}/p{discPor}LogCase{caseNo}.txt"
        totalLoadsFilePath = dataDir + f"/testMatrixTotal/curved/{loadCell}/p{discPor}LoadCase{caseNo}.txt"
        inertialLoadsFilePath = dataDir + f"/testMatrixInertial/{loadCell}/p{discPor}LoadCase{caseNo}.txt"
        
        self.totalLoadsObj = recordedData(totalLoadsFilePath, logFilePath)
        self.inertialLoadsObj = recordedData(inertialLoadsFilePath, logFilePath)
        self.totalLoadsObj.maxHarmonic = self.maxHarmonic
        self.inertialLoadsObj.maxHarmonic = self.maxHarmonic
        
        recMotorRPM = caseParams[int(caseNo),6]
        recPlungingRPM = recMotorRPM * 104 / 1000
        self.recSurgeFreq = recPlungingRPM / 60
        self.recSurgeFreqRed = 2 * np.pi * self.recSurgeFreq * 0.2 / self.totalLoadsObj.fstreamVel
        
        self.cycleQty = min(np.shape(self.totalLoadsObj.cycleDataClean)[0],np.shape(self.inertialLoadsObj.cycleDataClean)[0])
        print("Cycle Qty:", self.cycleQty)
        self.cycleSize = min(np.min(self.totalLoadsObj.cycleDataClean[:,1]),np.min(self.inertialLoadsObj.cycleDataClean[:,1]))
        
        self.resampledPhases = np.linspace(0,360,self.cycleSize,endpoint=False)
        
        self.totalLoadsObj.intpData, self.totalLoadsObj.intpCycleData = self.totalLoadsObj.resampleData(self.totalLoadsObj.dataExtdClean,self.totalLoadsObj.cycleDataClean,self.cycleQty,self.resampledPhases)
        self.inertialLoadsObj.intpData, self.inertialLoadsObj.intpCycleData = self.inertialLoadsObj.resampleData(self.inertialLoadsObj.dataExtdClean,self.inertialLoadsObj.cycleDataClean,self.cycleQty,self.resampledPhases)
        
        self.totalLoadsObj.filtIntpData, self.totalLoadsObj.filtIntpCycleData = self.totalLoadsObj.filterData(self.totalLoadsObj.intpData)
        self.inertialLoadsObj.filtIntpData, self.inertialLoadsObj.filtIntpCycleData = self.inertialLoadsObj.filterData(self.inertialLoadsObj.intpData)
        
        self.totalLoadsObj.filtData, self.totalLoadsObj.filtCycleData = self.totalLoadsObj.filterData(self.totalLoadsObj.dataExtdClean)
        self.inertialLoadsObj.filtData, self.inertialLoadsObj.filtCycleData = self.inertialLoadsObj.filterData(self.inertialLoadsObj.dataExtdClean)
        
        self.totalLoadsObj.intpFiltData, self.totalLoadsObj.intpFiltCycleData = self.totalLoadsObj.resampleData(self.totalLoadsObj.filtData,self.totalLoadsObj.filtCycleData,self.cycleQty,self.resampledPhases)
        self.inertialLoadsObj.intpFiltData, self.inertialLoadsObj.intpFiltCycleData = self.inertialLoadsObj.resampleData(self.inertialLoadsObj.filtData,self.inertialLoadsObj.filtCycleData,self.cycleQty,self.resampledPhases)
        
        self.aeroLoadsObj = self.initAeroObj()
        self.aeroLoadsObj.ndimIntpCycleAvgData, self.aeroLoadsObj.ndimIntpCycleStdData = self.aeroLoadsObj.cycleAveraging(self.aeroLoadsObj.ndimIntpCycleData)
        self.aeroLoadsObj.ndimDirectFiltCycleAvgData, self.aeroLoadsObj.ndimDirectFiltCycleStdData = self.aeroLoadsObj.cycleAveraging(self.aeroLoadsObj.ndimDirectFiltCycleData)
        self.aeroLoadsObj.ndimFiltIntpCycleAvgData, self.aeroLoadsObj.ndimFiltIntpCycleStdData = self.aeroLoadsObj.cycleAveraging(self.aeroLoadsObj.ndimFiltIntpCycleData)
        
        
    def initAeroObj(self):
        
        aeroLoadsObj = deepcopy(self.totalLoadsObj)
        del aeroLoadsObj.data
        del aeroLoadsObj.dataExtd
        del aeroLoadsObj.cycleData
        del aeroLoadsObj.dataExtdClean
        del aeroLoadsObj.cycleDataClean
        
        aeroLoadsObj.maxHarmonic = self.maxHarmonic
        
        aeroLoadsObj.intpData[:,2:] = self.totalLoadsObj.intpData[:,2:] - self.inertialLoadsObj.intpData[:,2:]
        aeroLoadsObj.directFiltData = aeroLoadsObj.filterData(aeroLoadsObj.intpData)[0]
        aeroLoadsObj.filtIntpData[:,2:] = self.totalLoadsObj.filtIntpData[:,2:] - self.inertialLoadsObj.filtIntpData[:,2:]
        aeroLoadsObj.intpFiltData[:,2:] = self.totalLoadsObj.intpFiltData[:,2:] - self.inertialLoadsObj.intpFiltData[:,2:]
        
        aeroLoadsObj.ndimIntpData = aeroLoadsObj.nonDimensionalize(aeroLoadsObj.intpData)
        aeroLoadsObj.ndimDirectFiltData = aeroLoadsObj.nonDimensionalize(aeroLoadsObj.directFiltData)
        aeroLoadsObj.ndimFiltIntpData = aeroLoadsObj.nonDimensionalize(aeroLoadsObj.filtIntpData)
        aeroLoadsObj.ndimIntpFiltData = aeroLoadsObj.nonDimensionalize(aeroLoadsObj.intpFiltData)
        
        aeroLoadsObj.ndimIntpCycleData = aeroLoadsObj.cycleSplit(aeroLoadsObj.ndimIntpData, isDataExtd=True)[1]
        aeroLoadsObj.ndimDirectFiltCycleData = aeroLoadsObj.cycleSplit(aeroLoadsObj.ndimDirectFiltData, isDataExtd=True)[1]
        aeroLoadsObj.ndimFiltIntpCycleData = aeroLoadsObj.cycleSplit(aeroLoadsObj.ndimFiltIntpData, isDataExtd=True)[1]
        aeroLoadsObj.ndimIntpFiltCycleData = aeroLoadsObj.cycleSplit(aeroLoadsObj.ndimIntpFiltData, isDataExtd=True)[1]
        
        return aeroLoadsObj
    
    def genAeroPlots(self):
        
        
        self.aeroLoadsObj.figIntpDFT, self.aeroLoadsObj.axIntpDFT = plt.subplots()
        dftPlot(self.aeroLoadsObj.ndimIntpData[:,2], 2e3, self.aeroLoadsObj.axIntpDFT, 20*self.aeroLoadsObj.surgeFreq, self.aeroLoadsObj.surgeFreq)
        
        self.aeroLoadsObj.figFiltDFT, self.aeroLoadsObj.axFiltDFT = plt.subplots()
        dftPlot(self.aeroLoadsObj.ndimDirectFiltData[:,2], 2e3, self.aeroLoadsObj.axFiltDFT, 20*self.aeroLoadsObj.surgeFreq, self.aeroLoadsObj.surgeFreq, color="tab:orange")
        
        # self.aeroLoadsObj.figFiltIntpDFT, self.aeroLoadsObj.axFiltIntpDFT = plt.subplots()
        # dftPlot(self.aeroLoadsObj.ndimFiltIntpData[:,2], 2e3, self.aeroLoadsObj.axFiltIntpDFT, 20*self.aeroLoadsObj.surgeFreq, self.aeroLoadsObj.surgeFreq, color="tab:green")
        
        self.aeroLoadsObj.figIntp, self.aeroLoadsObj.axIntp = plt.subplots()
        self.aeroLoadsObj.axIntp.plot(self.aeroLoadsObj.ndimIntpCycleAvgData[0], self.aeroLoadsObj.ndimIntpCycleAvgData[1])
        self.aeroLoadsObj.axIntp.fill_between(self.aeroLoadsObj.ndimIntpCycleAvgData[0],self.aeroLoadsObj.ndimIntpCycleAvgData[1] - self.aeroLoadsObj.ndimIntpCycleStdData[1], self.aeroLoadsObj.ndimIntpCycleAvgData[1] + self.aeroLoadsObj.ndimIntpCycleStdData[1],alpha=0.5)
        
        self.aeroLoadsObj.figFilt, self.aeroLoadsObj.axFilt = plt.subplots()
        self.aeroLoadsObj.axIntp.plot(self.aeroLoadsObj.ndimDirectFiltCycleAvgData[0], self.aeroLoadsObj.ndimDirectFiltCycleAvgData[1],color="tab:orange")
        self.aeroLoadsObj.axIntp.fill_between(self.aeroLoadsObj.ndimDirectFiltCycleAvgData[0],self.aeroLoadsObj.ndimDirectFiltCycleAvgData[1] - self.aeroLoadsObj.ndimDirectFiltCycleStdData[1], self.aeroLoadsObj.ndimDirectFiltCycleAvgData[1] + self.aeroLoadsObj.ndimDirectFiltCycleStdData[1],alpha=0.5, color="tab:orange")
        # self.aeroLoadsObj.axFilt.plot(self.aeroLoadsObj.ndimFiltIntpCycleAvgData[0], self.aeroLoadsObj.ndimFiltIntpCycleAvgData[1],color="tab:green")
        # self.aeroLoadsObj.axFilt.fill_between(self.aeroLoadsObj.ndimFiltIntpCycleAvgData[0],self.aeroLoadsObj.ndimFiltIntpCycleAvgData[1] - self.aeroLoadsObj.ndimFiltIntpCycleStdData[1], self.aeroLoadsObj.ndimFiltIntpCycleAvgData[1] + self.aeroLoadsObj.ndimFiltIntpCycleStdData[1],alpha=0.5, color="tab:green")
        
        
        
        
        
        
        
        
        