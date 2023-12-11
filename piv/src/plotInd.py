from matplotlib import pyplot as plt
import numpy as np
from ioUtilsNoStitch import loadData

class plotInd:
    
    def __init__(self, caseDataPath):
        
        self.data = loadData(caseDataPath)
        
        caseNo = int(self.data.case[-1])
        
        if caseNo < 3:
                
            staticFileName = caseDataPath.split("/")[-1].split(".")[0][:-1] + "6" + ".json"
            staticCaseDataPath = "/".join(caseDataPath.split("/")[:-1]) + "/" + staticFileName
            self.staticData = loadData(staticCaseDataPath)
            
        elif caseNo >= 3 and caseNo < 6:
            
            staticFileName = caseDataPath.split("/")[-1].split(".")[0][:-1] + "7" + ".json"
            staticCaseDataPath = "/".join(caseDataPath.split("/")[:-1]) + "/" + staticFileName
            self.staticData = loadData(staticCaseDataPath)
            
        elif caseNo == 6:
            
            dynamicFileName1 = caseDataPath.split("/")[-1].split(".")[0][:-1] + "0" + ".json"
            dynamicFileName2 = caseDataPath.split("/")[-1].split(".")[0][:-1] + "1" + ".json"
            dynamicFileName3 = caseDataPath.split("/")[-1].split(".")[0][:-1] + "2" + ".json"
            
            dynamicCaseDataPath1 = "/".join(caseDataPath.split("/")[:-1]) + "/" + dynamicFileName1
            dynamicCaseDataPath2 = "/".join(caseDataPath.split("/")[:-1]) + "/" + dynamicFileName2
            dynamicCaseDataPath3 = "/".join(caseDataPath.split("/")[:-1]) + "/" + dynamicFileName3
            
            self.dynamicData1 = loadData(dynamicCaseDataPath1)
            self.dynamicData2 = loadData(dynamicCaseDataPath2)
            self.dynamicData3 = loadData(dynamicCaseDataPath3)
            
        elif caseNo == 7:
            
            dynamicFileName1 = caseDataPath.split("/")[-1].split(".")[0][:-1] + "3" + ".json"
            dynamicFileName2 = caseDataPath.split("/")[-1].split(".")[0][:-1] + "4" + ".json"
            dynamicFileName3 = caseDataPath.split("/")[-1].split(".")[0][:-1] + "5" + ".json"
            
            dynamicCaseDataPath1 = "/".join(caseDataPath.split("/")[:-1]) + "/" + dynamicFileName1
            dynamicCaseDataPath2 = "/".join(caseDataPath.split("/")[:-1]) + "/" + dynamicFileName2
            dynamicCaseDataPath3 = "/".join(caseDataPath.split("/")[:-1]) + "/" + dynamicFileName3
            
            self.dynamicData1 = loadData(dynamicCaseDataPath1)
            self.dynamicData2 = loadData(dynamicCaseDataPath2)
            self.dynamicData3 = loadData(dynamicCaseDataPath3)
        
        #self.uid = f"k = {self.data.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{self.data.reducedSurgeAmplitude}"
        
    def phases(self, spanLocIndex, axsObj, plotStatic=False):
        
        #span location index between 0 and 20
        xData = self.data.ndimTimes
        
        yData = []
        
        for axialIndDistr in self.data.cycleAxialInd["intp200"]:
            yData.append(axialIndDistr[spanLocIndex])
            
        ndimSpanLoc = self.data.cycleAxialInd["spanLocs"][spanLocIndex]
        uid = f"r/R = {round(ndimSpanLoc,1)*2}"
        
        lines = []
        
        lines.append(axsObj.plot(xData, yData, marker="x", label=uid))
        
        if plotStatic:
            yStaticData = self.staticData.cycleAxialInd["intp200"][0][spanLocIndex]
            lines.append(axsObj.axhline(yStaticData, linestyle="--", color = lines[0][0].get_color()))
            
        return lines
        
    def distr(self, phaseNo, axsObj):
        
        ndimSpanLocs = np.array(self.data.cycleAxialInd["spanLocs"])*2
        yData = np.round(ndimSpanLocs,1)
        xData = self.data.cycleAxialInd["intp200"][0]
        xDynamicData1 = self.dynamicData1.cycleAxialInd["intp200"][phaseNo]
        xDynamicData2 = self.dynamicData2.cycleAxialInd["intp200"][phaseNo]
        xDynamicData3 = self.dynamicData3.cycleAxialInd["intp200"][phaseNo]
        
        label = f"k = {self.data.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{self.data.reducedSurgeAmplitude}"
        label1 = f"k = {self.dynamicData1.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{self.dynamicData1.reducedSurgeAmplitude}"
        label2 = f"k = {self.dynamicData2.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{self.dynamicData2.reducedSurgeAmplitude}"
        label3 = f"k = {self.dynamicData3.reducedSurgeFrequency}, " + r"$x_{amp}/D$ = " + f"{self.dynamicData3.reducedSurgeAmplitude}"
        
        axsObj.plot(xData, yData, marker="x", label=label, linestyle = "-", color="k")
        #axsObj.plot(xDynamicData1, yData, marker="x", label=label1)
        #axsObj.plot(xDynamicData2, yData, marker="x", label=label2)
        #axsObj.plot(xDynamicData3, yData, marker="x", label=label3)
        
        axsObj.grid()
        axsObj.legend()
        axsObj.set_xlabel(r"a [-]")
        axsObj.set_ylabel(r"r/R [-]")
            
    def surgeCycle6Span(self, plotDir):
        
        surgeCycleFig, surgeCycleAx = plt.subplots()
        
        centrePlot = self.phases(10,surgeCycleAx, plotStatic=True)
        radius0p2Plot = self.phases(12,surgeCycleAx, plotStatic=True)
        radius0p4Plot = self.phases(14,surgeCycleAx, plotStatic=True)
        radius0p6Plot = self.phases(16,surgeCycleAx, plotStatic=True)
        radius0p8Plot = self.phases(18,surgeCycleAx, plotStatic=True)
        radius1Plot = radius0p2Plot = self.phases(20,surgeCycleAx, plotStatic=True)
        
        surgeCycleAx.axvline(self.data.minLocTime, linestyle="dotted", color="red")
        surgeCycleAx.axvline(self.data.maxLocTime, linestyle="dotted", color="green")
        surgeCycleAx.axvline(self.data.meanLocTime1, linestyle="dotted", color="blue")
        surgeCycleAx.axvline(self.data.meanLocTime2, linestyle="dotted", color="blue")
        
        surgeCycleAx.grid()
        surgeCycleAx.legend()
        surgeCycleAx.set_xlabel(r"t/T [-]")
        surgeCycleAx.set_ylabel(r"a [-]")
        
        surgeCycleFig.savefig(plotDir + f"/surgeCycle6Span/{self.data.case}.png", dpi=300, bbox_inches="tight", pad_inches=0)
        
    def discPlane4Phase(self, plotDir):
        
        figMinLoc, axMinLoc = plt.subplots()
        figMaxLoc, axMaxLoc = plt.subplots()
        figMeanLoc1, axMeanLoc1 = plt.subplots()
        figMeanLoc2, axMeanLoc2 = plt.subplots()
        
        minLocPlot = self.distr(2, axMinLoc)
        maxLocPlot = self.distr(7, axMaxLoc)
        meanLoc1Plot = self.distr(4, axMeanLoc1)
        meanLoc2Plot = self.distr(9, axMeanLoc2)
        
        figMinLoc.savefig(plotDir + f"/discPlane4Phase/{self.data.case}MinLoc.png", dpi=300, bbox_inches="tight", pad_inches=0)
        figMaxLoc.savefig(plotDir + f"/discPlane4Phase/{self.data.case}MaxLoc.png", dpi=300, bbox_inches="tight", pad_inches=0)
        figMeanLoc1.savefig(plotDir + f"/discPlane4Phase/{self.data.case}MeanLoc1.png", dpi=300, bbox_inches="tight", pad_inches=0)
        figMeanLoc2.savefig(plotDir + f"/discPlane4Phase/{self.data.case}MeanLoc2.png", dpi=300, bbox_inches="tight", pad_inches=0)
        
        
        
        
        
        
        
            
        