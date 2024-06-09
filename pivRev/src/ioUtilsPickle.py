import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import CubicSpline
import pickle
from copy import deepcopy

from src.miscTools import genPolygonLineBndry
from src.miscTools import genVector
from src.miscTools import nrmlVector2D
from src.miscTools import lowPassFilt
from src.miscTools import movingAverage

class pickleData:
    
    def __init__(self, pickleDataDir, discPor, caseID):
        
        self.caseNo = int(caseID)
        self.caseDict = np.array([[1.0, 0.250, 3.0, 2.3906, 5.0e-2],
                             [2.0, 0.125, 3.0, 4.7830, 2.5e-2], 
                             [2.0, 0.250, 3.0, 4.7830, 5.0e-2], 
                             [2.0, 0.375, 2.0, 3.1875, 7.5e-2], 
                             [3.0, 0.250, 2.0, 4.7830, 5.0e-2],
                             [2.7, 0.375, 2.0, 4.3038, 7.5e-2],
                             [0.0, 0.000, 3.0, 0.0000, 0.0e-2],
                             [0.0, 0.000, 2.0, 0.0000, 0.0e-2]])
        
        with open(pickleDataDir + f"/p{discPor}Case{caseID}.pickle","rb") as handle:
        
            pickleObj = pickle.load(handle)
            
        for pickleAttrb in dir(pickleObj):
            
            if pickleAttrb[0] != "_":
                
                setattr(self,pickleAttrb,getattr(pickleObj,pickleAttrb))
                
        if self.caseNo < 6:
            
            period = 1/self.caseDict[self.caseNo,3]
            self.DT = 2*0.1*period
                
        self.genInterpolator()
        # self.getDiscVelRaw()
        # self.getDiscVelExtp()
        # self.pickleContainer = self.packData()
            
    def genInterpolator(self):
        
        x = self.gridX[0,:]
        y = self.gridY[:,0]
        
        for phaseObjID in self.phaseObjNames:
        
            phaseObj=getattr(self,phaseObjID)
            
            validCells = np.logical_not(phaseObj.gridMask)
            
            for gridAttrb in phaseObj.frameGridAttrbs[2:6]:
                
                # intpr = NearestNDInterpolator(np.column_stack((self.gridX[validCells],self.gridY[validCells])),getattr(phaseObj,gridAttrb)[validCells])
                intpr = RegularGridInterpolator((x,y),np.transpose(getattr(phaseObj,gridAttrb)))
                
                setattr(phaseObj,gridAttrb+"Intpr",intpr)
            
    def getDiscVelExtp(self, offset=0.20, attrb="Vx"):
        
        deltaX = self.gridX[0,1] - self.gridX[0,0]
        xArrSize = int(self.discDia/deltaX)
        
        # if self.caseNo < 6:
            
        #     offset = self.caseDict[self.caseNo,1]# Basing offset based on reduced amplitude
        #     print(offset)
        
        setattr(self,f"{attrb}cExtp",np.zeros((21,len(self.phaseObjNames))))
        setattr(self,f"xArr{attrb}cExtp",np.zeros((xArrSize,len(self.phaseObjNames))))
        
        polyfitObjs = np.empty((21,len(self.phaseObjNames)),dtype=np.ndarray)
        
        for i,phaseObjID in enumerate(self.phaseObjNames):
        
            phaseObj=getattr(self,phaseObjID)
        
            xEnd = phaseObj.discXc - offset * self.discDia
            xArr = np.linspace(xEnd-self.discDia,xEnd,xArrSize)
            getattr(self,f"xArr{attrb}cExtp")[:,i] =xArr
            
            yArr = np.linspace(phaseObj.discYc - 0.5*self.discDia,phaseObj.discYc + 0.5*self.discDia,21)
            discVelArr = np.zeros(len(yArr))
            
            for j,y in enumerate(yArr):
                
                line = np.zeros((len(xArr),2))
                line[:,0] = xArr
                line[:,1] = y
                
                vel = getattr(phaseObj,f"grid{attrb}Intpr")(line)
                polyfitObj= np.poly1d(np.polyfit(xArr,vel,deg=3))
                
                polyfitObjs[j,i] = polyfitObj
                discVelArr[j] = polyfitObj(phaseObj.discXc)
                
            getattr(self,f"{attrb}cExtp")[:,i] = discVelArr
            
        setattr(self,f"yArr{attrb}cExtp",yArr)
        
        return polyfitObjs
    
    def getDiscVelExtpV2(self, ustreamTrim=15, offset=0.25):
        
        polyfitObjs = np.empty((21,len(self.phaseObjNames)),dtype=np.ndarray)
        discVx = np.zeros((21,len(self.phaseObjNames)))
        
        for i,phaseObjID in enumerate(self.phaseObjNames):
            
            phaseObj=getattr(self,phaseObjID)
            
            xEnd = phaseObj.discXc - offset * self.discDia
            xEndIndex = np.argmin(np.abs(self.gridX[0,:] - xEnd))
            xArr = self.gridX[0,ustreamTrim:xEndIndex+1]
            
            yArr = np.linspace(phaseObj.discYc - 0.5*self.discDia,phaseObj.discYc + 0.5*self.discDia,21)
            
            for j,y in enumerate(yArr):
                
                line = np.zeros((len(xArr),2))
                line[:,0] = xArr
                line[:,1] = y
                
                vel = getattr(phaseObj,f"gridVxIntpr")(line)
                polyfitObj= np.poly1d(np.polyfit(xArr,vel,deg=3))
                
                polyfitObjs[j,i] = polyfitObj
                discVx[j,i] = polyfitObj(phaseObj.discXc)
                
        return polyfitObjs, discVx
    
    def getDiscVelRev(self,xcOffset=0):
        
        yOffsetArr = np.linspace(-0.5,0.5,21)
        discVxArr = np.zeros((21,len(self.phaseObjNames)))
        
        for i in range(21):
            
            velExtpObjs = self.getVxExtp(yOffsetArr[i],cutOff2=1)[1]
            
            for j in range(len(self.phaseObjNames)):
            
                discXc = getattr(self,self.phaseObjNames[j]).discXc
                discVxArr[i,j] = velExtpObjs[j](discXc + xcOffset*self.discDia)
                
        return discVxArr
    
    def getAvgDiscVel(self, annuliQty=100):
        
        yArr = np.linspace(-0.5,0.5,2*annuliQty + 1)
        temp = movingAverage(yArr)
        
        annuliLocs = temp[np.argwhere(temp >= 0)]
        annuliWidth = np.diff(yArr)[np.argwhere(temp >= 0)]
        annuliArea = 2*np.pi*annuliLocs*annuliWidth
        
        annuliAreaArr = np.zeros((annuliQty,len(self.phaseObjNames)))
        
        for j in range(len(self.phaseObjNames)):
            
            annuliAreaArr[:,j] = annuliArea.reshape(-1)
            
        Vx1Arr = np.zeros((annuliQty,len(self.phaseObjNames)))
        Vx2Arr = np.zeros((annuliQty,len(self.phaseObjNames)))
        
        for i,yLoc in enumerate(temp[:annuliQty]):
            
            velExtpObjs = self.getVxExtp(yLoc,cutOff2=1)[1]
            
            for j in range(len(self.phaseObjNames)):
                
                discXc = getattr(self,self.phaseObjNames[j]).discXc
                Vx1Arr[i,j] = velExtpObjs[j](discXc)
                
        for i,yLoc in enumerate(temp[annuliQty:]):
            
            velExtpObjs = self.getVxExtp(yLoc,cutOff2=1)[1]
            
            for j in range(len(self.phaseObjNames)):
                
                discXc = getattr(self,self.phaseObjNames[j]).discXc
                Vx2Arr[i,j] = velExtpObjs[j](discXc)
                
        Vx1Arr = np.flip(Vx1Arr,axis=0) #ensures that the first element corresponds to disc centre
        VxArr = 0.5 * (Vx1Arr + Vx2Arr)
        
        areaAvgNum = np.sum(VxArr * annuliAreaArr, axis=0)
        areaAvgDen = np.sum(annuliArea)
        
        avgDiscVel = areaAvgNum/areaAvgDen
        
        return avgDiscVel
    
    def getAvgDiscVelV2(self, annuliQty=100, offset=0.25):
        
        yArr = np.linspace(-0.5,0.5,2*annuliQty + 1)
        temp = movingAverage(yArr)
        
        annuliLocs = temp[np.argwhere(temp >= 0)]
        annuliWidth = np.diff(yArr)[np.argwhere(temp >= 0)]
        annuliArea = 2*np.pi*annuliLocs*annuliWidth
        
        annuliAreaArr = np.zeros((annuliQty,len(self.phaseObjNames)))
        
        for j in range(len(self.phaseObjNames)):
            
            annuliAreaArr[:,j] = annuliArea.reshape(-1)
            
        Vx1Arr = np.zeros((annuliQty,len(self.phaseObjNames)))
        Vx2Arr = np.zeros((annuliQty,len(self.phaseObjNames)))
        
        for i,yLoc in enumerate(temp[:annuliQty]):
            
            polyFitObjs = self.getVxExtpV2(yLoc, offset=offset)[1]
            
            for j in range(len(self.phaseObjNames)):
                
                discXc = getattr(self,self.phaseObjNames[j]).discXc
                Vx1Arr[i,j] = polyFitObjs[j](discXc)
                
        for i,yLoc in enumerate(temp[annuliQty:]):
            
            polyFitObjs = self.getVxExtpV2(yLoc, offset=offset)[1]
            
            for j in range(len(self.phaseObjNames)):
                
                discXc = getattr(self,self.phaseObjNames[j]).discXc
                Vx2Arr[i,j] = polyFitObjs[j](discXc)
                
        Vx1Arr = np.flip(Vx1Arr,axis=0) #ensures that the first element corresponds to disc centre
        VxArr = 0.5 * (Vx1Arr + Vx2Arr)
        
        areaAvgNum = np.sum(VxArr * annuliAreaArr, axis=0)
        areaAvgDen = np.sum(annuliArea)
        
        avgDiscVel = areaAvgNum/areaAvgDen
        
        return avgDiscVel
        
    def getVxExtp(self,yOffset=0,ustreamTrim=15,dstreamTrim=30,cutOff1=-0.25,cutOff2=0.5):
        
        '''Four points xStart, xEnd, xCutOff1, xCutOff2 are selected.
            xStart - Point close to inlet
            xEnd - Point close to outlet
            xCutOff1 - Point upstream of disc
            xCutOff2 - Point downstream of disc
            line - Line between xStart and xEnd
            line1 - Line between xStart and xCutOff1
            line2 - Line between xCutOff2 and xEnd
            line3 - Line between xCutOff1 and xCutOff2
            Vx at line 1 and line2 is directly obtained from the field
            Vx at line 3 is extrapolated from data at line 1 and line2'''
        
        xArr = self.gridX[0,ustreamTrim:-dstreamTrim+1]
        y = self.phase00.discYc + yOffset*self.discDia
        yArr = np.tile(y,len(xArr))
        linePts = np.column_stack((xArr, yArr))
        
        VxExtpObjs = np.empty((len(self.phaseObjNames)),dtype=np.ndarray)
        
        for i,phaseObjID in enumerate(self.phaseObjNames):
            
            phaseObj = getattr(self,phaseObjID)
            
            xCutOff1 = phaseObj.discXc + cutOff1 #Upstream cutoff point for Vx extrapolation
            xCutOff2 = phaseObj.discXc + cutOff2 #Downstream cutoff point for Vx extrapolation
            
            VxRaw = phaseObj.gridVxIntpr(linePts)
            VxFilt = lowPassFilt(VxRaw, len(VxRaw), 20) #Remove high freq noise to aid spline interpolation
            
            xCutOff1Indx = np.argmin(np.abs(linePts[:,0] - xCutOff1))
            xCutOff2Indx = np.argmin(np.abs(linePts[:,0] - xCutOff2))
            
            line1Pts = linePts[:xCutOff1Indx+1,:]
            line2Pts = linePts[xCutOff2Indx:,:]
            line12Pts = np.row_stack((line1Pts,line2Pts))
            
            VxFilt1 = VxFilt[:xCutOff1Indx+1]
            VxFilt2 = VxFilt[xCutOff2Indx:]
            VxFilt12 = np.concatenate((VxFilt1,VxFilt2))
            
            VxExtpObj = CubicSpline(line12Pts[:,0],VxFilt12)
            VxExtpObjs[i] = VxExtpObj
            
        return linePts, VxExtpObjs
    
    def getVxExtpV2(self,yOffset=0,ustreamTrim=15,offset=0.25):
        
        polyFitObjs = np.empty((len(self.phaseObjNames)),dtype=np.ndarray)
        
        for i,phaseObjID in enumerate(self.phaseObjNames):
            
            phaseObj = getattr(self,phaseObjID)
        
            xEnd = phaseObj.discXc - offset * self.discDia
            xEndIndex = np.argmin(np.abs(self.gridX[0,:] - xEnd))
            xArr = self.gridX[0,ustreamTrim:xEndIndex+1]
            
            y = self.phase00.discYc + yOffset*self.discDia
            yArr = np.tile(y,len(xArr))
            linePts = np.column_stack((xArr, yArr))
            
            vel = getattr(phaseObj,f"gridVxIntpr")(linePts)
            polyfitObj= np.poly1d(np.polyfit(xArr,vel,deg=3))
            polyFitObjs[i] = polyfitObj
            
        return linePts, polyFitObjs
            
            
            
    def getDiscVelRaw(self,offset=0.1):
        
        self.VxcRaw = np.zeros((21,len(self.phaseObjNames)))
        self.VycRaw = np.zeros((21,len(self.phaseObjNames)))
        self.VxuRaw = np.zeros((21,len(self.phaseObjNames)))
        self.VyuRaw = np.zeros((21,len(self.phaseObjNames)))
        
        for i,phaseObjID in enumerate(self.phaseObjNames):
        
            phaseObj=getattr(self,phaseObjID)
        
            xcLine = np.tile(phaseObj.discXc,21)
            xuLine = np.tile(phaseObj.discXc - offset*self.discDia,21)
            yLine = np.linspace(phaseObj.discYc - 0.5*self.discDia,phaseObj.discYc + 0.5*self.discDia,21)
            
            discCentreLine = np.column_stack((xcLine,yLine))
            discUstreamLine = np.column_stack((xuLine,yLine))
            
            self.VxcRaw[:,i] = phaseObj.gridVxIntpr(discCentreLine)
            self.VycRaw[:,i] = phaseObj.gridVyIntpr(discCentreLine)
            self.VxuRaw[:,i] = phaseObj.gridVxIntpr(discUstreamLine)
            self.VyuRaw[:,i] = phaseObj.gridVyIntpr(discUstreamLine)
            
        self.VcRaw = self.VxcRaw**2 + self.VycRaw**2
        self.VuRaw = self.VxuRaw**2 + self.VyuRaw**2
        
    def getCirc1D(self, pt1, pt2, VxExtpToggle=False):
        
        circ1D = np.zeros((len(self.phaseObjNames),1))
        lineBndryPtList = np.empty(len(self.phaseObjNames),dtype=np.ndarray)
         
        if VxExtpToggle == True:
             
            VxExtpObjs = self.getVxExtp(yOffset=pt1[1])[1]
         
        for i,phaseObjID in enumerate(self.phaseObjNames):
            
            phaseObj=getattr(self,phaseObjID)
            
            lineBndryPt1 = np.array([phaseObj.discXc + pt1[0]*self.discDia, phaseObj.discYc + pt1[1]*self.discDia])
            lineBndryPt2 = np.array([phaseObj.discXc + pt2[0]*self.discDia, phaseObj.discYc + pt2[1]*self.discDia])
            
            lineBndryPts = np.row_stack((lineBndryPt1,lineBndryPt2))
            
            bndryLineSeg, cntrlPtSeg = genPolygonLineBndry(lineBndryPts, lineDelta=self.gridX[0,1] - self.gridX[0,0])
            
            bndryVecSeg, distSeg = genVector(bndryLineSeg[0][:-1,:],bndryLineSeg[0][1:,:])
            
            if VxExtpToggle==False:
                VxSeg = phaseObj.gridVxIntpr(cntrlPtSeg[0])
            else:
                VxSeg = VxExtpObjs[i](cntrlPtSeg[0][:,0])
                
            VySeg = phaseObj.gridVyIntpr(cntrlPtSeg[0])
            
            circSeg = np.sum(np.column_stack((VxSeg, VySeg)) * bndryVecSeg * np.column_stack((distSeg, distSeg)))
            
            circ1D[i] = circSeg
            lineBndryPtList[i] = lineBndryPts
            
        return lineBndryPtList, circ1D
    
    def getWzFlux1D(self, pt1, pt2, VxExtpToggle=False):
        
        WzFlux1D = np.zeros((len(self.phaseObjNames),1))
        lineBndryPtList = np.empty(len(self.phaseObjNames),dtype=np.ndarray)
        
        if VxExtpToggle == True:
             
            VxExtpObjs = self.getVxExtp(yOffset=pt1[1])[1]
        
        for i,phaseObjID in enumerate(self.phaseObjNames):
            
            phaseObj=getattr(self,phaseObjID)
            
            lineBndryPt1 = np.array([phaseObj.discXc + pt1[0]*self.discDia, phaseObj.discYc + pt1[1]*self.discDia])
            lineBndryPt2 = np.array([phaseObj.discXc + pt2[0]*self.discDia, phaseObj.discYc + pt2[1]*self.discDia])
            
            lineBndryPts = np.row_stack((lineBndryPt1,lineBndryPt2))
            
            bndryLineSeg, cntrlPtSeg = genPolygonLineBndry(lineBndryPts, lineDelta=self.gridX[0,1] - self.gridX[0,0])
            
            bndryVecSeg, distSeg = nrmlVector2D(bndryLineSeg[0][:-1,:],bndryLineSeg[0][1:,:])[1:]
            
            if VxExtpToggle==False:
                VxSeg = phaseObj.gridVxIntpr(cntrlPtSeg[0])
            else:
                VxSeg = VxExtpObjs[i](cntrlPtSeg[0][:,0])
                
            VySeg = phaseObj.gridVyIntpr(cntrlPtSeg[0])
            
            WzSeg = phaseObj.gridWzIntpr(cntrlPtSeg[0])
            
            WzFluxSeg = np.sum(np.column_stack((WzSeg,WzSeg)) * np.column_stack((VxSeg, VySeg)) * bndryVecSeg * np.column_stack((distSeg, distSeg)))
            
            WzFlux1D[i] = WzFluxSeg
            lineBndryPtList[i] = lineBndryPts
            
        return lineBndryPtList, WzFlux1D
     
    def getCircWzFlux2D(self, xBounds=np.array([-0.25,1]), yBounds=np.array([-0.8,0])):
         
        pt1 = [xBounds[0], yBounds[0]]
        pt2 = [xBounds[1], yBounds[0]]
        pt3 = [xBounds[1], yBounds[1]]
        pt4 = [xBounds[0], yBounds[1]]
        
        if np.abs(yBounds[0]) <= 0.5:
            line1Pts, line1Circ = self.getCirc1D(pt1, pt2, VxExtpToggle=True)
            line1WzFlux = self.getWzFlux1D(pt1, pt2, VxExtpToggle=True)[1]
        else:
            line1Pts, line1Circ = self.getCirc1D(pt1, pt2)
            line1WzFlux = self.getWzFlux1D(pt1, pt2)[1]
            
        line2Pts, line2Circ = self.getCirc1D(pt2, pt3)
        line2WzFlux = self.getWzFlux1D(pt2, pt3)[1]
        
        if np.abs(yBounds[0]) <= 0.5:
            line3Pts, line3Circ = self.getCirc1D(pt3, pt4, VxExtpToggle=True)
            line3WzFlux = self.getWzFlux1D(pt3, pt4, VxExtpToggle=True)[1]
        else:
            line3Pts, line3Circ = self.getCirc1D(pt3, pt4)
            line3WzFlux = self.getWzFlux1D(pt3, pt4)[1]
            
        line4Pts, line4Circ = self.getCirc1D(pt4, pt1)
        line4WzFlux = self.getWzFlux1D(pt4, pt1)[1]
        
        lineSegPts = np.column_stack((line1Pts, line2Pts, line3Pts, line4Pts))
        lineSegCirc = np.column_stack((line1Circ, line2Circ, line3Circ, line4Circ))
        lineSegWzFlux = np.column_stack((line1WzFlux, line2WzFlux, line3WzFlux, line4WzFlux))
         
        return lineSegPts, lineSegCirc, lineSegWzFlux
    
    def circWzParamStudy(self, xMin=-0.3, xMax=1.5, yMin=-0.8, yMax=0.8, deltaX=0.1, deltaY=0.05):
        
        xU = xMin
        if xU < 0.3:
            xDArr = np.arange(0.3,xMax+deltaX,deltaX)
        else:
            xDArr = np.arange(xU+0.3,xMax+deltaX,deltaX)
        yLArr = np.arange(-0.45,deltaY,deltaY)
        yUArr = np.arange(0.45,-deltaY,-deltaY)
        
        lowerData = np.empty((len(xDArr),len(yLArr)),dtype=np.ndarray)
        upperData = np.empty((len(xDArr),len(yUArr)),dtype=np.ndarray)
        
        for i in range(len(xDArr)):
            for j in range(len(yLArr)):
                
                lowerData[i,j] = cntrlVolContainer(self,[xU,xDArr[i]],[yMin,yLArr[j]])
                
        for i in range(len(xDArr)):
            for j in range(len(yUArr)):
                
                upperData[i,j] = cntrlVolContainer(self,[xU,xDArr[i]],[yUArr[j],yMax])
        
        return lowerData, upperData
    
    def getCntrlVolArrData(self, cntrlVolArr, attrb="gamma", plotPhaseID="00"):
        
        attrbData = np.zeros((cntrlVolArr.shape))
        
        for i in range(cntrlVolArr.shape[0]):
            for j in range(cntrlVolArr.shape[1]):
                
                attrbData[i,j] = getattr(cntrlVolArr[i,j],attrb)[int(plotPhaseID)]
                
        return attrbData
                
    def getCirc2D(self, xBounds=np.array([-1,1]), yBounds=np.array([-0.85,0.85])):
        
        self.circ = np.zeros((len(self.phaseObjNames),4))
        rectBndryPtList = np.empty(len(self.phaseObjNames),dtype=np.ndarray)
        
        for i,phaseObjID in enumerate(self.phaseObjNames):
            
            phaseObj=getattr(self,phaseObjID)
            
            rectBndryPt1 = np.array([phaseObj.discXc + xBounds[0]*self.discDia, phaseObj.discYc + yBounds[0]*self.discDia])
            rectBndryPt2 = np.array([phaseObj.discXc + xBounds[1]*self.discDia, phaseObj.discYc + yBounds[0]*self.discDia])
            rectBndryPt3 = np.array([phaseObj.discXc + xBounds[1]*self.discDia, phaseObj.discYc + yBounds[1]*self.discDia])
            rectBndryPt4 = np.array([phaseObj.discXc + xBounds[0]*self.discDia, phaseObj.discYc + yBounds[1]*self.discDia])
            
            rectBndryPts = np.row_stack((rectBndryPt1,rectBndryPt2,rectBndryPt3,rectBndryPt4))
            
            bndryLineSegs, cntrlPtSegs = genPolygonLineBndry(rectBndryPts, lineDelta=self.gridX[0,1] - self.gridX[0,0])
            
            circSegs = np.zeros(len(bndryLineSegs))
            
            for j in range(len(bndryLineSegs)):
                
                bndryVecSeg, distSeg = genVector(bndryLineSegs[j][:-1,:],bndryLineSegs[j][1:,:])
                VxSeg = phaseObj.gridVxIntpr(cntrlPtSegs[j])
                VySeg = phaseObj.gridVyIntpr(cntrlPtSegs[j])
                
                circSegs[j] = np.sum(np.column_stack((VxSeg, VySeg)) * bndryVecSeg * np.column_stack((distSeg, distSeg)))
                
            self.circ[i,:] = circSegs
            rectBndryPtList[i] = rectBndryPts
            
        return rectBndryPtList
        
    def maskFrame(self, phaseObj):
        
        maskPhaseObj = deepcopy(phaseObj)
        
        for gridAttrb in phaseObj.frameGridAttrbs[2:-1]:
            
            setattr(maskPhaseObj,gridAttrb,np.ma.masked_array(getattr(phaseObj,gridAttrb),phaseObj.gridMask))
            
        return maskPhaseObj
    
    def packData(self):
        
        return caseContainer(self)
    
class cntrlVolContainer:
    
    def __init__(self, pickleObj, xLimits, yLimits):
        
        self.xLimits = xLimits
        self.yLimits = yLimits
        self.caseNo = pickleObj.caseNo
        
        gammaSegs, omegaFluxSegs = pickleObj.getCircWzFlux2D(self.xLimits,self.yLimits)[1:]
        # gamma = np.sum(gamma,axis=1)
        # omegaFlux = np.sum(omegaFlux,axis=1)
        
        if  self.caseNo < 6:
            # gammaTemp = np.append(gamma[10:],np.append(gamma[2:],gamma[[2,3]]))
            # DGamma_DTTemp = (np.roll(gammaTemp,-1) - (np.roll(gammaTemp,1)))/pickleObj.DT
            gammaTemp = np.row_stack((gammaSegs[10:,:],gammaSegs[2:,:],gammaSegs[[2,3],:]))
            DGamma_DTTemp = np.gradient(gammaTemp,0.1,axis=0,edge_order=2)
            
            self.attrs = ["gammaSegs", "DGamma_DTSegs", "omegaFluxSegs"]
            
            # self.gamma = np.append(gamma[2:],gamma[2])
            # self.DGamma_DT = DGamma_DTTemp[2:-1]
            # self.omegaFlux = np.append(omegaFlux[2:],omegaFlux[2])
            
            self.gammaSegs = np.row_stack((gammaSegs[2:,:],gammaSegs[2,:]))
            self.DGamma_DTSegs = DGamma_DTTemp[2:-1,:]
            self.omegaFluxSegs = np.row_stack((omegaFluxSegs[2:,:],omegaFluxSegs[2,:]))
            
        else:
            
            self.attrs = ["gammaSegs", "omegaFluxSegs"]
            
            self.gammaSegs = gammaSegs
            self.omegaFluxSegs = omegaFluxSegs
        
class caseContainer:
    
    def __init__(self, pickleObj):
        
        self.uid = pickleObj.uid
        self.discDia = pickleObj.discDia
        
        self.gridX = pickleObj.gridX
        self.gridY = pickleObj.gridY
        
        self.phaseObjNames = pickleObj.phaseObjNames
        
        for phaseName in pickleObj.phaseObjNames:
        
            setattr(self,phaseName,phaseContainer(getattr(pickleObj,phaseName),pickleObj))
            
        self.getDiscVelRev = pickleObj.getDiscVelRev.__func__
        print("Generating control volume arrays...")
        self.lowerCVArr, self.upperCVArr = pickleObj.circWzParamStudy(deltaX=0.01,deltaY=0.01)
        print("Control volume arrays generated!")
        self.getCntrlVolArrData = pickleObj.getCntrlVolArrData.__func__
            
        # self.discVxcExtp = pickleObj.VxcExtp
        # self.discVxuRaw = pickleObj.VxuRaw
                    
class phaseContainer:
    
    def __init__(self, phaseObj, pickleObj):
        
        maskedField = pickleObj.maskFrame(phaseObj)
        
        self.rhoInf = phaseObj.rhoInf
        self.VInf = phaseObj.VInf
        self.discXc = phaseObj.discXc
        self.discYc = phaseObj.discYc
        self.frameGridAttrbs = phaseObj.frameGridAttrbs
        
        for gridAttrb in maskedField.frameGridAttrbs[2:-1]:
            
            setattr(self,gridAttrb,getattr(maskedField,gridAttrb))
        
            
            
            