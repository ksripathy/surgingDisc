import numpy as np
import glob
from pivFields import planarPIVField
from dataclasses import dataclass
from dataclasses_json import dataclass_json

#Use ioUtils file for data in tecPlotStitch with only one combined siby by side piv frame data

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
        p45DiscLoc = np.array([[-1.151, -24.425, -32.129, -22.294, 2.291, 34.745, 61.297, 70.148, 58.511, 30.319, -1.642, -24.917],
                               [6.717, -4.265, -7.707, -0.659, 13.273, 29.991, 42.284, 45.890, 38.023, 23.435, 6.553, -4.101],
                               [-3.282, -26.228, -32.293, -23.770, 2.291, 35.564, 65.067, 70.640, 60.642, 30.647, -3.773, -26.884],
                               [-8.445, -41.554, -54.667, -44.013, -9.265, 40.234, 83.997, 98.585, 83.178, 41.710, -10.084, -43.521],
                               [-0.086, -23.852, -29.097, -20.738, 5.487, 38.431, 68.098, 73.671, 64.001, 33.678, -0.578, -23.361],
                               [-8.937, -46.307, -51.880, -44.668, -9.921, 41.382, 87.767, 102.030, 85.144, 41.710, -10.576, -45.816]])
                              
        p70DiscLoc = np.array([[-0.742, -24.672, -31.884, -22.213, 2.536, 34.498, 61.214, 69.901, 57.936, 29.908, -1.397, -24.836],
                               [6.962, -4.184, -7.134, -0.414, 13.026, 30.072, 41.710, 45.807, 37.776, 23.680, 6.634, -4.012],
                               [-3.364, -27.131, -31.720, -23.852, 2.372, 35.153, 64.165, 69.246, 60.559, 29.581, -3.692, -27.294],
                               [-5.207, -38.687, -52.827, -38.687, -4.395, 43.225, 88.082, 103.20, 86.781, 43.225, -6.020, -39.337],
                               [3.082, -20.322, -25.685, -16.909, 8.932, 41.600, 71.342, 77.192, 66.953, 36.724, 2.431, -20.322],
                               [-6.680, -43.328, -50.876, -40.312, -5.857, 44.850, 90.682, 105.150, 87.269, 44.688, -7.645, -42.263]])
        
        staticDiscLoc = 21.058
        
        self.uid = tecFile.split("/")[-1].split(".")[0]
        self.discPor = float(self.uid[1:3])/100
        caseNo = int(self.uid[7])
        self.redFreq = caseDict[caseNo, 0]
        self.redAmp = caseDict[caseNo, 1]
        self.fstreamVelc = caseDict[caseNo, 2]
        
            
        for file in fileList:
                        
            phaseID = int(file.split(".")[0].split("_")[1])
            setattr(self, f"phase{phaseID}", planarPIVField(file))
            
        #self.axialIndDiscMid = self.axialInd()
        
        #assigning axial loication of disc
        if self.discPor == 0.45 and self.redFreq != 0.0:
            
            discLocs = p45DiscLoc[caseNo,:]
            
        elif self.discPor == 0.70 and self.redFreq != 0.0:
            
            discLocs = p70DiscLoc[caseNo,:]
            
        elif self.redFreq == 0.0:
            
            discLocs = np.array([staticDiscLoc])
            
        #assigning spanwise location of disc. Since P70 cases 3,4,5 are performed after camera reinstallation the spanwise location is different due to revised calibration method
        if self.discPor == 0.70 and (caseNo == 3 or caseNo == 4 or caseNo == 5):
             
            self.axialIndDiscMid = self.axialIndv2(7.0, discLocs)
             
        else:
             
            self.axialIndDiscMid = self.axialIndv2(72.0, discLocs)
            
        #Compute event times for dynamic cases only
        if self.redFreq != 0.0:
            self.eventTimes(discLocs)    
           
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
        arch.axialInductionDiscCentre = self.axialIndDiscMid
        
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
    axialInductionDiscCentre: np.float64 = np.zeros(phases)
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