import numpy as np
import glob
from src.pivFields import planarPIVField
from dataclasses import dataclass, field
from dataclasses_json import dataclass_json
from typing import Dict
from scipy.io import savemat
from src.miscTools import nearestValueIndex

class caseData:
    
    def __init__(self, tecFile):
        
        fileList = glob.glob(tecFile)
        self.phaseQty = len(fileList)
        
        #Lookup table for case metadata
        caseDict = np.array([[1.0, 0.250, 3.0],
                             [2.0, 0.125, 3.0], 
                             [2.0, 0.250, 3.0], 
                             [2.0, 0.375, 2.0], 
                             [3.0, 0.250, 2.0],
                             [2.7, 0.375, 2.0],
                             [0.0, 0.000, 3.0],
                             [0.0, 0.000, 2.0]])
        
        #Lookup table for disc location
        p45DiscLoc = np.array([[3.855, -19.480, -27.204, -17.344, 7.470, 40.008, 66.630, 75.350, 63.672, 35.407, 3.362, -19.973],
                               [12.236, 0.733, -2.554, 4.512, 18.316, 35.078, 47.403, 51.018, 43.130, 28.669, 11.743, 0.897],
                               [1.883, -21.617, -27.368, -18.659, 7.470, 40.665, 70.081, 75.997, 65.808, 35.900, 1.390, -22.110],
                               [-3.075, -36.559, -49.008, -38.735, -3.568, 45.403, 89.116, 103.740, 88.294, 46.882, -4.390, -38.571],
                               [3.082, -20.322, -25.685, -16.909, 8.932, 41.600, 71.342, 77.192, 66.953, 36.724, 2.431, -20.322],
                               [-3.897, -42.022, -47.166, -39.064, -4.390, 46.553, 93.060, 107.190, 90.266, 46.882, -5.047, -41.693]])
                              
        p70DiscLoc = np.array([[-0.742, -24.672, -31.884, -22.213, 2.536, 34.498, 61.214, 69.901, 57.936, 29.908, -1.397, -24.836],
                               [6.962, -4.184, -7.134, -0.414, 13.026, 30.072, 41.710, 45.807, 37.776, 23.680, 6.634, -4.012],
                               [-3.364, -27.131, -31.720, -23.852, 2.372, 35.153, 64.165, 69.246, 60.559, 29.581, -3.692, -27.294],
                               [-5.207, -38.687, -52.827, -38.687, -4.395, 43.225, 88.082, 103.20, 86.781, 43.225, -6.020, -39.337],
                               [3.082, -20.322, -25.685, -16.909, 8.932, 41.600, 71.342, 77.192, 66.953, 36.724, 2.431, -20.322],
                               [-6.680, -43.328, -50.876, -40.312, -5.857, 44.850, 90.682, 105.150, 87.269, 44.688, -7.645, -42.263]])
        
        p45DiscCentreLoc = np.array([5.932, 5.932, 5.932, 8.069, 8.069, 8.069])
        p70DiscCentreLoc = np.array([8.069, 8.069, 8.069, 6.369, 6.369, 6.369])
        
        p45StaticDiscLoc = 25.683
        p70StaticDiscLoc = 26.012
        
        staticDiscCentreLoc = 8.069
        
        self.uid = tecFile.split("/")[-1].split(".")[0]
        self.case=self.uid[:8]
        self.discPor = float(self.uid[1:3])/100
        self.caseNo = int(self.uid[7])
        self.redFreq = caseDict[self.caseNo, 0]
        self.redAmp = caseDict[self.caseNo, 1]
        self.fstreamVelc = caseDict[self.caseNo, 2]
        self.discDia = 200
        
        #assigning axial loication of disc
        if self.discPor == 0.45 and self.redFreq != 0.0:
            
            discAxialLocs = p45DiscLoc[self.caseNo,:]
            discCentreLoc = p45DiscCentreLoc[self.caseNo]
            
        elif self.discPor == 0.70 and self.redFreq != 0.0:
            
            discAxialLocs = p70DiscLoc[self.caseNo,:]
            discCentreLoc = p70DiscCentreLoc[self.caseNo]
            
        elif self.discPor == 0.45 and self.redFreq == 0.0:
            
            discAxialLocs = np.array([p45StaticDiscLoc])
            discCentreLoc = staticDiscCentreLoc
            
        elif self.discPor == 0.70 and self.redFreq == 0.0:
            
            discAxialLocs = np.array([p70StaticDiscLoc])
            discCentreLoc = staticDiscCentreLoc
                
        self.matDict = {"case" : self.uid[:8], "k" : self.redFreq, "amp" : self.redAmp, 
                       "Vinf" : self.fstreamVelc, "rcIndex" : 0, "xcIndex" : np.zeros(len(discAxialLocs))}
        
        #Compute event times for dynamic cases only
        if self.redFreq != 0.0:
            self.eventTimes(discAxialLocs)
            self.matDict["phaseTimes"] = self.cumElapsedTimes        
            
        for i, file in enumerate(fileList):
                        
            phaseID = int(file.split(".")[0].split("_")[1])
            temp = planarPIVField(file)
            temp.combine2Frames()
            
            setattr(temp.combinedFrames,"phaseVinf",np.abs(temp.combinedFrames.recFstreamVelc()))
            
            frameAxialLocs = temp.frameAxialLocs
            frameSpanLocs = temp.frameSpanLocs
            
            self.matDict["deltaX"] = (frameAxialLocs[1] - frameAxialLocs[0]) * 1e-3
            self.matDict["deltaY"] = (frameSpanLocs[0] - frameSpanLocs[1]) * 1e-3
            self.matDict["rcIndex"] = nearestValueIndex(frameSpanLocs, discCentreLoc) + 1
            self.matDict["xcIndex"][i] = nearestValueIndex(frameAxialLocs, discAxialLocs[i]) + 1 
            self.matDict |= {f"phase{phaseID+1}" : {}}
            self.matDict[f"phase{phaseID+1}"] |= {"X" : temp.combinedFrames.gridPosX*1e-3}
            self.matDict[f"phase{phaseID+1}"] |= {"Y" : temp.combinedFrames.gridPosY*1e-3}
            self.matDict[f"phase{phaseID+1}"] |= {"Vx" : temp.combinedFrames.gridVelX}
            self.matDict[f"phase{phaseID+1}"] |= {"Vr" : temp.combinedFrames.gridVelY}
            self.matDict[f"phase{phaseID+1}"] |= {"Vmag" : temp.combinedFrames.gridVel}
            self.matDict[f"phase{phaseID+1}"] |= {"Vinf" : temp.combinedFrames.phaseVinf}
            
            setattr(self, f"phase{phaseID}", temp)
            
        #self.axialIndDistr = self.axialInd()
            
        #self.cycleAxialInd = self.cycleAxialIndDistr(discCentreLoc, discAxialLocs)
            
    def saveAsMat(self, fileDir):
        
        savemat(fileDir + "/" + self.uid[:8] + ".mat", self.matDict)
           
    def axialInd(self):

        axialIndDiscMid = np.zeros(self.phaseQty)
        
        for i in range(self.phaseQty):
            
            pivObj = getattr(self, f"phase{i}")
            axialIndDiscMid[i] = pivObj.discAxialInd()[1]
            
        return axialIndDiscMid
    
    def axialIndv2(self, spanLoc, discLocs):
        
        axialIndDiscMid = np.zeros(self.phaseQty)
        
        for i in range(self.phaseQty):
            
            pivObj = getattr(self, f"phase{i}")
            axialIndDiscMid[i] = pivObj.discAxialIndv2(spanLoc, discLocs[i])
            
        return axialIndDiscMid
    
    def cycleAxialIndDistr(self, discCentreLoc, discAxialLocs, distrSize=21):
        
        cycleAxialIndDistrIntp0  = np.empty(self.phaseQty, dtype=np.ndarray)
        cycleAxialIndDistrIntp200  = np.empty(self.phaseQty, dtype=np.ndarray)
        
        for phaseNo, axialLoc in enumerate(discAxialLocs):
            
            pivObj = getattr(self, f"phase{phaseNo}")
            pivObj.combine2Frames(overlapSmoothing=True)
            
            cycleAxialIndDistrIntp0[phaseNo] = pivObj.discAxialIndDistrV2(1, discCentreLoc, axialLoc, 0,distrSize)
            cycleAxialIndDistrIntp200[phaseNo] = pivObj.discAxialIndDistrV2(1, discCentreLoc, axialLoc, 200, distrSize)
            
        return {"spanLocs": (np.linspace(discCentreLoc - self.discDia*0.5, discCentreLoc + self.discDia*0.5, distrSize) - discCentreLoc)/ self.discDia,
                "intp0": cycleAxialIndDistrIntp0,
                "intp200": cycleAxialIndDistrIntp200}
    
    def eventTimes(self, discLocs):
        
        meanPos = 0.5 * (np.min(discLocs) + np.max(discLocs))
        discLocsScaled = discLocs - meanPos
        dispAmp = 0.5 * (np.max(discLocsScaled) - np.min(discLocsScaled))
        #Time taken from sin motion equation. Here mean position time corresponds to 0
        #Min position time corresponds to -0.25 and max position time corresponds to 0.25 
        computedTimes = np.arcsin(np.clip(discLocsScaled/dispAmp,-1,1))/(2*np.pi)
        
        #Elapsed times between phases
        elapsedTimes = np.abs(computedTimes[1:] - np.roll(computedTimes,1)[1:])
        
        #Cumulative elapsedTimes
        self.cumElapsedTimes = np.zeros(12)
        self.cumElapsedTimes[1:] = np.cumsum(elapsedTimes)
        
        self.minLocTime = self.cumElapsedTimes[2]
        self.maxLocTime = self.cumElapsedTimes[7]
        self.meanLocTime1 = self.cumElapsedTimes[4] + np.abs(computedTimes[4])
        self.meanLocTime2 = self.cumElapsedTimes[9] + np.abs(computedTimes[9])
    
    def save(self, fileDir):
        
        arch = caseDataArch()
        
        arch.case = self.uid[:8]
        arch.discPorosity = self.discPor
        arch.reducedSurgeFrequency = self.redFreq
        arch.reducedSurgeAmplitude = self.redAmp
        arch.freestreamVelocity = self.fstreamVelc
        arch.phases = self.phaseQty
        arch.cycleAxialInd = self.cycleAxialInd
        
        if self.redFreq != 0:
        
            arch.ndimTimes = self.cumElapsedTimes
            arch.minLocTime = self.minLocTime
            arch.maxLocTime = self.maxLocTime
            arch.meanLocTime1 = self.meanLocTime1
            arch.meanLocTime2 = self.meanLocTime2
        
        arch.saveAsJson(fileDir)           
            
#Serialize case data into json archive
@dataclass_json
@dataclass
class caseDataArch:
    
    case: str = "Case data archive not initialized"
    discPorosity: float = 0.0
    reducedSurgeFrequency: float = 0.0
    reducedSurgeAmplitude: float = 0.0
    freestreamVelocity: float = 0.0
    phases: int = 12
    spanLocs: np.float64 = np.zeros(phases)
    cycleAxialInd: Dict = field(default_factory=dict) #Initialization syntax for mutable objects like dict and list
    ndimTimes: np.float64 = np.zeros(phases)
    minLocTime: float = 0.0
    maxLocTime: float = 0.0
    meanLocTime1: float = 0.0
    meanLocTime2: float = 0.0
    
    def saveAsJson(self, fileDir):
        
        with open(f"{fileDir}/results_{self.case}.json", 'w') as jsonFile:
            
            jsonFile.write(self.to_json(indent=4))


def loadData(fileLoc):
    
    with open(fileLoc, 'r') as jsonFile:
        
        input = caseDataArch.from_json(jsonFile.read())
        
    return input