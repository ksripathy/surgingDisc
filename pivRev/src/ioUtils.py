import numpy as np
import glob
from src.pivFields import planarPIVField

class caseData:
    
    def __init__(self, dataDir, discPor, caseID):
        
        caseNo = int(caseID)
        self.uid = f"p{discPor}Case{caseID}"
        print(f"Processing {self.uid}...\n")
        
        #Lookup table for case metadata
        caseDict = np.array([[1.0, 0.250, 3.0, 2.3906, 5.0e-2],
                             [2.0, 0.125, 3.0, 4.7830, 2.5e-2], 
                             [2.0, 0.250, 3.0, 4.7830, 5.0e-2], 
                             [2.0, 0.375, 2.0, 3.1875, 7.5e-2], 
                             [3.0, 0.250, 2.0, 4.7830, 5.0e-2],
                             [2.7, 0.375, 2.0, 4.3038, 7.5e-2],
                             [0.0, 0.000, 3.0, 0.0000, 0.0e-2],
                             [0.0, 0.000, 2.0, 0.0000, 0.0e-2]])
        
        self.redFreq = caseDict[caseNo,0]
        self.redAmp = caseDict[caseNo,1]
        self.redSurgeVel = self.redFreq * self.redAmp
        self.freq = caseDict[caseNo,3]
        self.amp = caseDict[caseNo,4]
        
        p45DiscAxialLocs = np.array([[3.855, -19.480, -27.204, -17.344, 7.470, 40.008, 66.630, 75.350, 63.672, 35.407, 3.362, -19.973],
                               [12.236, 0.733, -2.554, 4.512, 18.316, 35.078, 47.403, 51.018, 43.130, 28.669, 11.743, 0.897],
                               [1.883, -21.617, -27.368, -18.659, 7.470, 40.665, 70.081, 75.997, 65.808, 35.900, 1.390, -22.110],
                               [-3.075, -36.559, -49.008, -38.735, -3.568, 45.403, 89.116, 103.740, 88.294, 46.882, -4.390, -38.571],
                               [3.082, -20.322, -25.685, -16.909, 8.932, 41.600, 71.342, 77.192, 66.953, 36.724, 2.431, -20.322],
                               [-3.897, -42.022, -47.166, -39.064, -4.390, 46.553, 93.060, 107.190, 90.266, 46.882, -5.047, -41.693],
                               [25.683, None, None, None, None, None, None, None, None, None, None, None],
                               [25.683, None, None, None, None, None, None, None, None, None, None, None]])

        p70DiscAxialLocs = np.array([[-0.742, -24.672, -31.884, -22.213, 2.536, 34.498, 61.214, 69.901, 57.936, 29.908, -1.397, -24.836],
                               [6.962, -4.184, -7.134, -0.414, 13.026, 30.072, 41.710, 45.807, 37.776, 23.680, 6.634, -4.012],
                               [-3.364, -27.131, -31.720, -23.852, 2.372, 35.153, 64.165, 69.246, 60.559, 29.581, -3.692, -27.294],
                               [-5.207, -38.687, -52.827, -38.687, -4.395, 43.225, 88.082, 103.20, 86.781, 43.225, -6.020, -39.337],
                               [3.082, -20.322, -25.685, -16.909, 8.932, 41.600, 71.342, 77.192, 66.953, 36.724, 2.431, -20.322],
                               [-6.680, -43.328, -50.876, -40.312, -5.857, 44.850, 90.682, 105.150, 87.269, 44.688, -7.645, -42.263],
                               [26.012, None, None, None, None, None, None, None, None, None, None, None],
                               [26.012, None, None, None, None, None, None, None, None, None, None, None]])

        p45DiscCentreLocs = np.array([5.932, 5.932, 5.932, 8.069, 8.069, 8.069, 8.069, 8.069])
        p70DiscCentreLocs = np.array([8.069, 8.069, 8.069, 6.369, 6.369, 6.369, 8.069, 8.069])
        
        self.caseDataPaths = sorted(glob.glob(dataDir + f"/p{discPor}MeanCase{caseID}*"))
        self.caseLogPaths = sorted(glob.glob(dataDir + f"/p{discPor}LogCase{caseID}*"))
        
        if discPor == "45":
    
            discAxialLocs = p45DiscAxialLocs
            discCentreLocs = p45DiscCentreLocs
    
        elif discPor == "70":
            
            discAxialLocs = p70DiscAxialLocs
            discCentreLocs = p70DiscCentreLocs
            
        self.discXcArr = discAxialLocs[caseNo,:]
        self.discYc = discCentreLocs[caseNo]
        
        self.phaseObjNames = []
        
        for i, (phaseDataPath, phaseLogPath) in enumerate(zip(self.caseDataPaths,self.caseLogPaths)):
            
            phaseObjName = "phase" + phaseDataPath.split(".")[0][-2:]
            self.phaseObjNames.append(phaseObjName)
            print(f"Processing {phaseObjName}...")
            phaseObj = planarPIVField(phaseDataPath,phaseLogPath,self.discXcArr[i],self.discYc)
            print(f"Finished Processing {phaseObjName}!\n")
            
            setattr(self,phaseObjName,phaseObj)
            
        # self.genIndVel()
        # self.pickleContainer = self.packData()
        # self.pickleContainer = self.packFilledData()
        self.pickleContainer = self.packRawData()
        
        print(f"Finished Processing {self.uid}!\n")
        
    def genIndVel(self):
        
        indVelName = ["Vx", "Vy", "V"]
        indVelLoc = ["u","c"]
        indVelSrc = ["Raw", "Intp"]
        
        #Store induced velocities
        for name in indVelName:
            
            for loc in indVelLoc:
                
                for src in indVelSrc:
                    
                    indVel = np.zeros((21,len(self.phaseObjNames)))
                    
                    for i,phaseName in enumerate(self.phaseObjNames):
                        
                        phaseObj = getattr(self,phaseName)
                        
                        indVel[:,i] = getattr(phaseObj,"ind"+name+loc+src)
                        
                    setattr(self,"ind"+name+loc+src+"Arr",indVel)
                    
        #Store locations
        self.yLine = self.phase00.discCentreLine[:,1]
        
        self.xcArr = np.zeros((len(self.phaseObjNames)))
        self.xuArr = np.zeros((len(self.phaseObjNames)))
        
        for i,phaseName in enumerate(self.phaseObjNames):
            
            phaseObj = getattr(self,phaseName)
            
            self.xcArr[i] = phaseObj.discCentreLine[0,0]
            self.xuArr[i] = phaseObj.discUstreamLine[0,0]
            
    def packData(self):
        
        return caseContainerMinimal(self)
    
    def packFilledData(self):
        
        return filledCaseContainer(self)
    
    def packRawData(self):
        
        return caseContainerRaw(self)
            
class caseContainer:
    
    def __init__(self, ioUtilsObj):
        
        self.uid = ioUtilsObj.uid
        self.caseDataPaths = ioUtilsObj.caseDataPaths
        self.caseDataPaths = ioUtilsObj.caseLogPaths
        
        self.gridX = ioUtilsObj.phase00.ndimField.gridX
        self.gridY = ioUtilsObj.phase00.ndimField.gridY
        
        self.phaseObjNames = ioUtilsObj.phaseObjNames
        
        for phaseName in ioUtilsObj.phaseObjNames:
        
            setattr(self,phaseName,phaseContainer(getattr(ioUtilsObj,phaseName)))
            
        self.yLine = ioUtilsObj.yLine
        self.xcArr = ioUtilsObj.xcArr
        self.xuArr = ioUtilsObj.xuArr
        self.indVelLinear = indVelContainer(ioUtilsObj, "Raw")
        self.indVelRBF = indVelContainer(ioUtilsObj, "Intp")
        
class caseContainerRaw:
    
    def __init__(self, ioUtilsObj):
        
        self.uid = ioUtilsObj.uid
        
        self.gridX = ioUtilsObj.phase00.field.gridX
        self.gridY = ioUtilsObj.phase00.field.gridY
        
        self.discDia = ioUtilsObj.phase00.field.discDia
        self.phaseObjNames = ioUtilsObj.phaseObjNames
        
        for phaseName in ioUtilsObj.phaseObjNames:
        
            setattr(self,phaseName,phaseContainerRaw(getattr(ioUtilsObj,phaseName)))
        
class caseContainerMinimal:
    
    def __init__(self, ioUtilsObj):
        
        self.uid = ioUtilsObj.uid
        self.caseDataPaths = ioUtilsObj.caseDataPaths
        self.caseDataPaths = ioUtilsObj.caseLogPaths
        
        self.gridX = ioUtilsObj.phase00.ndimField.gridX
        self.gridY = ioUtilsObj.phase00.ndimField.gridY
        
        self.discDia = ioUtilsObj.phase00.ndimField.discDia
        
        self.phaseObjNames = ioUtilsObj.phaseObjNames
        
        for phaseName in ioUtilsObj.phaseObjNames:
        
            setattr(self,phaseName,phaseContainerMinimal(getattr(ioUtilsObj,phaseName)))
            
class filledCaseContainer:
    
    def __init__(self, ioUtilsObj):
        
        self.uid = ioUtilsObj.uid
        self.caseDataPaths = ioUtilsObj.caseDataPaths
        self.caseDataPaths = ioUtilsObj.caseLogPaths
        
        self.gridX = ioUtilsObj.phase00.ndimFieldFilled.gridX
        self.gridY = ioUtilsObj.phase00.ndimFieldFilled.gridY
        
        self.discDia = ioUtilsObj.phase00.ndimFieldFilled.discDia
        
        self.phaseObjNames = ioUtilsObj.phaseObjNames
        
        for phaseName in ioUtilsObj.phaseObjNames:
        
            setattr(self,phaseName,filledPhaseContainer(getattr(ioUtilsObj,phaseName)))
            
class phaseContainer:
    
    def __init__(self, phaseObj):
        
        maskedField = phaseObj.ndimField.maskFrame()
        
        self.rhoInf = phaseObj.rhoInf
        self.VInf = phaseObj.VInf
        self.frameGridAttrbs = phaseObj.frameGridAttrbs
        
        for gridAttrb in maskedField.frameGridAttrbs[2:-1]:
            
            setattr(self,gridAttrb,getattr(maskedField,gridAttrb))
            
            
        # self.gridVxIntpr = phaseObj.ndimField.gridVxIntpr
        # self.gridVyIntpr = phaseObj.ndimField.gridVyIntpr
        # self.gridWzIntpr = phaseObj.ndimField.gridWzIntpr
        # self.gridVxRBFIntpr = phaseObj.ndimField.gridVxRBFIntpr
        # self.gridVyRBFIntpr = phaseObj.ndimField.gridVyRBFIntpr
        
class phaseContainerMinimal:
    
    def __init__(self, phaseObj):
        
        frameObj = phaseObj.ndimField
        
        self.rhoInf = phaseObj.rhoInf
        self.VInf = phaseObj.VInf
        self.discXc = frameObj.discXc
        self.discYc = frameObj.discYc
        self.frameGridAttrbs = frameObj.frameGridAttrbs
        
        for gridAttrb in frameObj.frameGridAttrbs[2:]:
            
            setattr(self,gridAttrb,getattr(frameObj,gridAttrb))
            
        self.gridVxIntpr = phaseObj.ndimField.gridVxIntpr
        self.gridVyIntpr = phaseObj.ndimField.gridVyIntpr
        self.gridWzIntpr = phaseObj.ndimField.gridWzIntpr
        
class phaseContainerRaw:

    def __init__(self, phaseObj):
        
        frameObj = phaseObj.field
        
        self.rhoInf = phaseObj.rhoInf
        self.VInf = phaseObj.VInf
        self.discXc = frameObj.discXc
        self.discYc = frameObj.discYc
        self.frameGridAttrbs = frameObj.frameGridAttrbs
        
        for gridAttrb in frameObj.frameGridAttrbs[2:]:
            
            setattr(self,gridAttrb,getattr(frameObj,gridAttrb))
        
class filledPhaseContainer:
    
    def __init__(self, phaseObj):
        
        frameObj = phaseObj.ndimFieldFilled
        
        self.rhoInf = phaseObj.rhoInf
        self.VInf = phaseObj.VInf
        self.discXc = frameObj.discXc
        self.discYc = frameObj.discYc
        self.frameGridAttrbs = frameObj.frameGridAttrbs
        
        for gridAttrb in frameObj.frameGridAttrbs[2:]:
            
            setattr(self,gridAttrb,getattr(frameObj,gridAttrb))
            
        self.discMask = frameObj.discMask
            
class indVelContainer:
    
    def __init__(self, ioUtilsObj, srcType):
        
        indVelName = ["Vx", "Vy", "V"]
        indVelLoc = ["u","c"]
        
        #Store induced velocities
        for name in indVelName:
            
            for loc in indVelLoc:
            
                setattr(self,name+loc,getattr(ioUtilsObj,"ind"+name+loc+srcType+"Arr"))
            
        
            
            
                    
            
            
        