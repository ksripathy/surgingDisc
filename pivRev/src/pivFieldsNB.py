import numpy as np
from copy import deepcopy
from scipy.interpolate import interpn
from scipy.interpolate import griddata
from scipy.stats import linregress 
import matplotlib.pyplot as plt
import io
from scipy.interpolate import CubicSpline
from src.pivFrames import pivFrames2D
# from piv.src.miscTools import nearestValueIndex

class planarPIVField():
    
    def __init__(self, meanFilePath, logFilePath, discXc, discYc):
        
        self.framesCols, self.framesRows, self.framesTecData = self.readTecFile(meanFilePath, logFilePath)
        
    def readTecFile(self,meanFilePath,logFilePath):
        
        logData = np.loadtxt(logFilePath)
        
        self.rhoInf = logData[10]
        self.VInf = logData[12]
        
        with open(meanFilePath, "r") as tecFile:
            tecFileLines = tecFile.read().split("\n")
            
        frameStartIndxs = []
        
        for i,line in enumerate(tecFileLines):
            
            if "ZONE" in line.split():
                frameStartIndxs.append(i)
                
        framesData = np.empty(len(frameStartIndxs),dtype=np.ndarray)
        framesCols = np.empty(len(frameStartIndxs),dtype=int)
        framesRows = np.empty(len(frameStartIndxs),dtype=int)
                
        for i, frameStartIndx in enumerate(frameStartIndxs):
            
            frameCols, frameRows = self.getFrameSize(tecFileLines,frameStartIndx)
            frameSize = frameCols * frameRows
            
            frameFileLines = tecFileLines[frameStartIndx+2:frameStartIndx+2+frameSize]
            
            framesData[i] = np.loadtxt(io.StringIO("\n".join(frameFileLines)))
            framesCols[i] = frameCols
            framesRows[i] = frameRows
            
        return framesCols, framesRows, framesData   
                
    def getFrameSize(self,tecFileLines,frameStartIndx):
        
        temp = tecFileLines[frameStartIndx].split(",")
        
        frameWidth = int(temp[1].split("=")[1])
        frameHeight = int(temp[2].split("=")[1])
        
        return frameWidth, frameHeight
    
    def combineFrames(self, frameObj1, frameObj2, ovlpToggle=False, expFactor=12):
        
        #initialising combined frame as frameObj2
        combinedFrameObj = deepcopy(frameObj2)
        
        frameObj1.validCells = np.logical_not(frameObj1.gridMask)
        frameObj2.validCells = np.logical_not(frameObj2.gridMask)
        
        #Assigning frameObj1 data in non-ovlap region
        for gridAttrb in combinedFrameObj.frameGridAttrbs[2:]:
            
            getattr(combinedFrameObj,gridAttrb)[np.logical_and(frameObj1.validCells,np.logical_not(frameObj2.validCells))] = getattr(frameObj1,gridAttrb)[np.logical_and(frameObj1.validCells,np.logical_not(frameObj2.validCells))]
            
        if ovlpToggle:
            
            #Assigning data in overlap region
            frameOvlpIndxs = np.argwhere(np.logical_and(frameObj1.validCells,frameObj2.validCells))
            xOvlpMin = min(frameOvlpIndxs[:,1])
            xOvlpMax = max(frameOvlpIndxs[:,1])
            xOvlpNorm = (combinedFrameObj.gridX[frameOvlpIndxs[:,0],frameOvlpIndxs[:,1]] - xOvlpMin)/(xOvlpMax - xOvlpMin)
            
            frameObj1.ovlpWeights = 1/(1 + np.exp(expFactor * (xOvlpNorm - 0.5)))
            frameObj2.ovlpWeights = 1/(1 + np.exp(-expFactor * (xOvlpNorm - 0.5)))
            
            for gridAttrb in combinedFrameObj.frameGridAttrbs[2:-1]:
                
                getattr(combinedFrameObj,gridAttrb)[frameOvlpIndxs[:,0],frameOvlpIndxs[:,1]] = frameObj1.ovlpWeights*getattr(frameObj1,gridAttrb)[frameOvlpIndxs[:,0],frameOvlpIndxs[:,1]] + frameObj2.ovlpWeights * getattr(frameObj2,gridAttrb)[frameOvlpIndxs[:,0],frameOvlpIndxs[:,1]]

        combinedFrameObj.flattenData()
        
        return combinedFrameObj
    
    def getDiscMask(self,frameObj):
        
        discXc = frameObj.discXc
        discYc = frameObj.discYc
        discDia = frameObj.discDia
        
        point1 = (discXc + 0.35 * discDia, discYc + 0.75 * discDia)
        point2 = (discXc - 0.25 * discDia, discYc + 0.75 * discDia)
        point3 = (discXc - 0.25 * discDia, discYc - 0.75 * discDia)
        point4 = (discXc + 0.35 * discDia, discYc - 0.75 * discDia)
        
        polygonCheck = frameObj.getPolygonIndxs([point1, point2, point3, point4])
        
        discMaskCheck = np.logical_and(polygonCheck,frameObj.gridMask)
            
        return polygonCheck, discMaskCheck
    
    def getRBFMask(self,frameObj):
        
        discXc = frameObj.discXc
        discYc = frameObj.discYc
        discDia = frameObj.discDia
        
        point1 = (discXc + 0.35 * discDia, discYc + 0.85 * discDia)
        point2 = (discXc - 0.35 * discDia, discYc + 0.85 * discDia)
        point3 = (discXc - 0.35 * discDia, discYc - 0.85 * discDia)
        point4 = (discXc + 0.35 * discDia, discYc - 0.85 * discDia)
        
        polygonCheck = frameObj.getPolygonIndxs([point1, point2, point3, point4])
        
        return polygonCheck
    
    def getRBFMask2(self,frameObj,ustrOffset=0,dstrOffset=0,ustrmPolygonWidth=0.75,dstrmPolygonWidth=0.75):
        
        discXc = frameObj.discXc
        discYc = frameObj.discYc
        discDia = frameObj.discDia
        
        diskMask = self.getDiscMask(frameObj)[1]
        diskMaskIndxs = np.argwhere(diskMask)
        
        ustrmX1 = frameObj.gridX[0,min(diskMaskIndxs[:,1])] - ustrOffset * discDia
        ustrmX2 = ustrmX1 - ustrmPolygonWidth * discDia 
        dstrmX1 = frameObj.gridX[0,max(diskMaskIndxs[:,1])] + dstrOffset * discDia
        dstrmX2 = dstrmX1 + dstrmPolygonWidth * discDia
        
        minY = discYc - 0.85 * discDia
        maxY = discYc + 0.95 * discDia
        
        ustrmPolygonCheck = frameObj.getPolygonIndxs([(ustrmX1,maxY),(ustrmX2,maxY),(ustrmX2,minY),(ustrmX1,minY)])
        dstrmPolygonCheck = frameObj.getPolygonIndxs([(dstrmX2,maxY),(dstrmX1,maxY),(dstrmX1,minY),(dstrmX2,minY)])
        polygonSumCheck = frameObj.getPolygonIndxs([(dstrmX2,maxY),(ustrmX2,maxY),(ustrmX2,minY),(dstrmX2,minY)])
        polygonDiffCheck = frameObj.getPolygonIndxs([(dstrmX1,maxY),(ustrmX1,maxY),(ustrmX1,minY),(dstrmX1,minY)])
        
        return ustrmPolygonCheck, dstrmPolygonCheck, polygonSumCheck, polygonDiffCheck
    
    def getRBFMask3(self,frameObj,ustrOffset=0,dstrOffset=0,ustrmPolygonWidth=0.75,dstrmPolygonWidth=0.75):
        
        discXc = frameObj.discXc
        discYc = frameObj.discYc
        discDia = frameObj.discDia
        
        diskMask = self.getDiscMask(frameObj)[1]
        diskMaskIndxs = np.argwhere(diskMask)
        
        ustrmX1 = discXc - ustrOffset * discDia
        ustrmX2 = ustrmX1 - ustrmPolygonWidth * discDia 
        dstrmX1 = discXc + dstrOffset * discDia
        dstrmX2 = dstrmX1 + dstrmPolygonWidth * discDia
        
        minY = discYc - 0.85 * discDia
        maxY = discYc + 0.90 * discDia
        
        ustrmPolygonCheck = frameObj.getPolygonIndxs([(ustrmX1,maxY),(ustrmX2,maxY),(ustrmX2,minY),(ustrmX1,minY)])
        dstrmPolygonCheck = frameObj.getPolygonIndxs([(dstrmX2,maxY),(dstrmX1,maxY),(dstrmX1,minY),(dstrmX2,minY)])
        polygonSumCheck = frameObj.getPolygonIndxs([(dstrmX2,maxY),(ustrmX2,maxY),(ustrmX2,minY),(dstrmX2,minY)])
        polygonDiffCheck = frameObj.getPolygonIndxs([(dstrmX1,maxY),(ustrmX1,maxY),(ustrmX1,minY),(dstrmX1,minY)])
        
        return ustrmPolygonCheck, dstrmPolygonCheck, polygonSumCheck, polygonDiffCheck
    
    def getIndVelRaw(self, line):
        
        indVxRaw = self.ndimField.gridVxIntpr(line)
        indVyRaw = self.ndimField.gridVyIntpr(line)
        indVRaw = np.sqrt(indVxRaw**2 + indVyRaw**2)
        
        return indVxRaw, indVyRaw, indVRaw
    
    def getIndVelIntp(self,line):
        
        indVxIntp = self.ndimField.gridVxRBFIntpr(line)
        indVyIntp = self.ndimField.gridVyRBFIntpr(line)
        indVIntp = np.sqrt(indVxIntp**2 + indVyIntp**2)
        
        return indVxIntp, indVyIntp, indVIntp
    
    def getDiscVel(self, frameObj, offset=0.1, velAttrb="gridVx"):
        
        xEnd = frameObj.discXc - offset * frameObj.discDia
        xArr = np.linspace(xEnd-frameObj.discDia,xEnd,np.round(frameObj.discDia/frameObj.deltaX))
        
        yArr = np.linspace(frameObj.discYc - 0.5*frameObj.discDia,frameObj.discYc + 0.5*frameObj.discDia,21)
        discVelArr = np.zeros(len(xArr))
        
        for i,y in enumerate(yArr):
            
            line = np.zeros((2,len(xArr)))
            line[:,0] = xArr
            line[:,1] = y
            
            vel = frameObj.getDataLine(line,gridAttrb=velAttrb)
            splineObj = CubicSpline(xArr,line,bc_type="natural",extrapolate=True)
            
            discVelArr[i] = splineObj(frameObj.discXc)
            
        return discVelArr, xArr
    
    def getIndVelExtp(self, offset=0.1):
        
        indVxExtp = self.getDiscVel(self.ndimField,offset=offset,velAttrb="gridVx")
        indVyExtp = self.getDiscVel(self.ndimField,offset=offset,velAttrb="gridVy")
        indVExtp = np.sqrt(indVxExtp**2 + indVyExtp**2)
        
        return indVxExtp, indVyExtp, indVExtp
    
    def fillDiscMask(self, frameObj):
        
        filledFrameObj = deepcopy(frameObj)
        discMask = self.getDiscMask(frameObj)[1]
        
        frameGridAttrbs = deepcopy(filledFrameObj.frameGridAttrbs)
        
        del frameGridAttrbs[4:6]
        
        xMask = filledFrameObj.gridX[discMask]
        yMask = filledFrameObj.gridY[discMask]
        maskPts = np.column_stack((xMask,yMask))
        
        for gridAttrb in frameGridAttrbs[2:-1]:
            
            getattr(filledFrameObj,gridAttrb)[discMask] = getattr(frameObj,gridAttrb+"Intpr")(maskPts)
            
        filledFrameObj.gridV[discMask] = np.sqrt(filledFrameObj.gridVx[discMask]**2 + filledFrameObj.gridVy[discMask]**2)
        filledFrameObj.Wz = filledFrameObj.computeWz()
        filledFrameObj.gridMask[discMask] = False
        filledFrameObj.discMask = discMask
        
        filledFrameObj.flattenData()
        
        return filledFrameObj 
    
    
        
        
            
        
        
        
        
        
        
        

    
    
        
        
        
        
                   
            
            

            