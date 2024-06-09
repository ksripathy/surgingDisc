import numpy as np
from copy import deepcopy
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import CloughTocher2DInterpolator
from scipy.interpolate import RBFInterpolator
from scipy.stats import linregress
from numpy.ma import masked_array
from shapely import Polygon
from shapely import Point

class pivFrames2D:
    
    def __init__(self, colQty, rowQty, tecData, discXc, discYc, tecStdData=None):
        
        self.frameAttrbs = ["x", "y", "Vx", "Vy", "V", "Wz", "stdVx", "stdVy", "mask"]
        self.frameGridAttrbs = ["gridX", "gridY", "gridVx", "gridVy", "gridV", "gridWz", "gridStdVx", "gridStdVy", "gridMask"]
        
        self.x = tecData[:,0]
        self.y = tecData[:,1]
        self.Vx = tecData[:,2]
        self.Vy = tecData[:,3]
        self.V = np.sqrt(self.Vx**2 + self.Vy**2)
        self.Wz = tecData[:,9]
        self.stdVx = tecData[:,-3]
        self.stdVy = tecData[:,-2]
        self.mask = np.logical_not(tecData[:,-1].astype(np.bool_))
        
        self.gridData(colQty,rowQty)
        
        self.deltaX = np.round(self.gridX[0,1] - self.gridX[0,0],3)
        self.deltaY = np.round(self.gridY[0,0] - self.gridY[1,0],3)
        
        self.discXc = discXc
        self.discYc = discYc
        self.discDia = 200
        
    def gridData(self,colQty,rowQty):
        
        for attrb,gridAttrb in zip(self.frameAttrbs,self.frameGridAttrbs):
            
            attrbData = getattr(self,attrb)
            attrbGridData = attrbData.reshape((rowQty,colQty))
            setattr(self,gridAttrb,attrbGridData)
            
    def flattenData(self):
        
        for attrb, gridAttrb in zip(self.frameAttrbs,self.frameGridAttrbs):
            
            setattr(self,attrb,getattr(self,gridAttrb).flatten())
            
    def shiftFrame(self, srcFrame, colShift=0, rowShift=0):
        
        destFrame = deepcopy(srcFrame)
        
        destFrame.gridX = np.roll(srcFrame.gridX, colShift, axis=1)
        destFrame.gridY = np.roll(srcFrame.gridY, rowShift, axis=0)
        
        #Correcting boundary artifacting from numpy roll
        
        if colShift > 0:
            
            for i in np.flip(np.arange(colShift)):
                
                destFrame.gridX[:,i] = destFrame.gridX[:,i+1] - destFrame.deltaX
                
        elif colShift < 0:
            
            for i in np.arange(colShift,0):
                
                destFrame.gridX[:,i] = destFrame.gridX[:,i-1] + destFrame.deltaX
                
        if rowShift > 0:
            
            for i in np.flip(np.arange(rowShift)):
                
                destFrame.gridY[i,:] = destFrame.gridY[i+1,:] + destFrame.deltaY
                
        elif rowShift < 0:
            
            for i in np.arange(rowShift,0):
                
                destFrame.gridY[i,:] = destFrame.gridY[i-1,:] - destFrame.deltaY
                
        destFrame.flattenData()
                
        return destFrame
    
    def maskFrame(self):
        
        maskFrameObj = deepcopy(self)
        
        for gridAttrb in self.frameGridAttrbs[2:-1]:
            
            setattr(maskFrameObj,gridAttrb,masked_array(getattr(self,gridAttrb),self.mask))
            
        return maskFrameObj
    
    def genInterpolator(self):
        
        validCells = np.logical_not(self.gridMask)
        
        # for gridAttrb in self.frameGridAttrbs[2:5]:
        for gridAttrb in ["gridVx", "gridVy", "gridStdVx", "gridStdVy"]:
            
            intpr = LinearNDInterpolator(np.column_stack((self.gridX[validCells],self.gridY[validCells])),getattr(self,gridAttrb)[validCells])
            
            setattr(self,gridAttrb+"Intpr",intpr)
    
    def genRBFInterpolator(self, polygonMask, skipDims=1):
        
        maskIndxs = np.argwhere(polygonMask)
        
        rowMin = min(maskIndxs[:,0])
        rowMax = max(maskIndxs[:,0]) + 1
        colMin = min(maskIndxs[:,1])
        colMax = max(maskIndxs[:,1]) + 1
        
        validCells = np.logical_not(self.gridMask[rowMin:rowMax:skipDims,colMin:colMax:skipDims])
        
        print(f"Data length = {len(self.gridMask[rowMin:rowMax:skipDims,colMin:colMax:skipDims][validCells])}")
        
        for gridAttrb in self.frameGridAttrbs[2:4]:
    
            intpr = RBFInterpolator(np.column_stack((self.gridX[rowMin:rowMax:skipDims,colMin:colMax:skipDims][validCells],
                                                     self.gridY[rowMin:rowMax:skipDims,colMin:colMax:skipDims][validCells])),
                                    getattr(self,gridAttrb)[rowMin:rowMax:skipDims,colMin:colMax:skipDims][validCells])
            print("RBF Interpolant construction finished for " + gridAttrb)
            
            setattr(self,gridAttrb+"RBFIntpr",intpr)
            
    def genRBFInterpolator2(self, ustrmPolygonMask, dstrmPolygonMask, skipDims=1, rbfKernel="thin_plate_spline"):
        
        ustrmMaskIndxs = np.argwhere(ustrmPolygonMask)
        dstrmMaskIndxs = np.argwhere(dstrmPolygonMask)
        
        ustrmRowMin = min(ustrmMaskIndxs[:,0])
        ustrmRowMax = max(ustrmMaskIndxs[:,0]) + 1
        ustrmColMin = min(ustrmMaskIndxs[:,1])
        ustrmColMax = max(ustrmMaskIndxs[:,1]) + 1
        
        dstrmRowMin = min(dstrmMaskIndxs[:,0])
        dstrmRowMax = max(dstrmMaskIndxs[:,0]) + 1
        dstrmColMin = min(dstrmMaskIndxs[:,1])
        dstrmColMax = max(dstrmMaskIndxs[:,1]) + 1
        
        validUstrmCells = np.logical_not(self.gridMask[ustrmRowMin:ustrmRowMax:skipDims,ustrmColMin:ustrmColMax:skipDims])
        validDstrmCells = np.logical_not(self.gridMask[dstrmRowMin:dstrmRowMax:skipDims,dstrmColMin:dstrmColMax:skipDims])
        
        ustmX = self.gridX[ustrmRowMin:ustrmRowMax:skipDims,ustrmColMin:ustrmColMax:skipDims][validUstrmCells]
        ustmY = self.gridY[ustrmRowMin:ustrmRowMax:skipDims,ustrmColMin:ustrmColMax:skipDims][validUstrmCells]
        
        dstmX = self.gridX[dstrmRowMin:dstrmRowMax:skipDims,dstrmColMin:dstrmColMax:skipDims][validDstrmCells]
        dstmY = self.gridY[dstrmRowMin:dstrmRowMax:skipDims,dstrmColMin:dstrmColMax:skipDims][validDstrmCells]
        
        x = np.concatenate((ustmX,dstmX))
        y = np.concatenate((ustmY,dstmY))
        
        # print(f"X-length = {len(x)}, Y-Length = {len(y)}")
        
        for gridAttrb in self.frameGridAttrbs[2:4]:
            
            ustmZ = getattr(self,gridAttrb)[ustrmRowMin:ustrmRowMax:skipDims,ustrmColMin:ustrmColMax:skipDims][validUstrmCells]
            dstmZ = getattr(self,gridAttrb)[dstrmRowMin:dstrmRowMax:skipDims,dstrmColMin:dstrmColMax:skipDims][validDstrmCells]
            z = np.concatenate((ustmZ,dstmZ))
            
            print(f"X-length = {len(x)}, Y-Length = {len(y)}, Z-Length = {len(z)}")
    
            intpr = RBFInterpolator(np.column_stack((x,y)),z,kernel=rbfKernel)
            print("RBF Interpolant construction finished for " + gridAttrb)
            
            setattr(self,gridAttrb+"RBFIntpr",intpr)        
            
    def computeWz(self):
        
        dVy_dx = np.gradient(self.gridVy,self.gridX[0,:]*1e-3,axis=1)
        dVx_dy = np.gradient(self.gridVx,self.gridY[:,0]*1e-3,axis=0)
        
        return dVy_dx - dVx_dy
                
    def computeWz2(self):
        
        '''Computation of vorticity without edge artifacts arising from masks'''
        
        dv_dx = np.zeros(self.gridX.shape)
        du_dy = np.zeros(self.gridX.shape)
        
        #Scanning through each row of the field for dv_dx
        for i in range(self.gridX.shape[0]):
            
            validIndxs = np.argwhere(np.logical_not(self.gridMask[i,:])).reshape(-1)
            
            #discontinuity locs
            disCtys = np.argwhere(np.diff(validIndxs) != 1).reshape(-1)
            
            if len(validIndxs) > 3 and len(disCtys) == 0:
                
                dv_dx[i,validIndxs] = np.gradient(self.gridVy[i,validIndxs],self.gridX[i,validIndxs],edge_order=2)
            
            #If there are discontinuities between the valid cells of row    
            elif len(disCtys) > 0:
                
                #Compute voticity for upstream valid cells
                if disCtys[0]+1 - validIndxs[0] > 3:
                    
                    jSlice = np.arange(validIndxs[0],disCtys[0]+1)
                    dv_dx[i,jSlice] = np.gradient(self.gridVy[i,jSlice],self.gridX[i,jSlice],edge_order=2)
                    
                #Compute vorticity for downstream valid cells
                if validIndxs[-1]+1 - disCtys[-1] > 3:
                    
                    jSlice = np.arange(disCtys[-1],validIndxs[-1]+1)
                    dv_dx[i,jSlice] = np.gradient(self.gridVy[i,jSlice],self.gridX[i,jSlice],edge_order=2)
                    
                #Compute vorticity valid cells other than theabove two
                if len(disCtys) > 1:
                    
                    for disCtyIndx in range(1,len(disCtys)-1):
                        
                        validIndx1 = validIndxs[disCtys[disCtyIndx]+1]
                        validIndx2 = validIndxs[disCtys[disCtyIndx+1]]
                        
                        if validIndx2+1 - validIndx1 > 3:
                            
                            jSlice = np.arange(validIndx1,validIndx2+1)
                            dv_dx[i,jSlice] = np.gradient(self.gridVy[i,jSlice],self.gridX[i,jSlice],edge_order=2)
                
        #Scanning through each row of the field for du_dy
        for j in range(self.gridX.shape[1]):
            
            validIndxs = np.argwhere(np.logical_not(self.gridMask[:,j])).reshape(-1)
            
            #discontinuity locs
            disCtys = np.argwhere(np.diff(validIndxs) != 1).reshape(-1)
            
            if len(validIndxs) > 3 and len(disCtys) == 0:
                
                du_dy[validIndxs,j] = np.gradient(self.gridVx[validIndxs,j],self.gridY[validIndxs,j],edge_order=2)
                
            
            elif len(disCtys) > 0:
                
                if disCtys[0]+1 - validIndxs[0] > 3:
                    
                    iSlice = np.arange(validIndxs[0],disCtys[0]+1)
                    du_dy[iSlice,j] = np.gradient(self.gridVx[iSlice,j],self.gridY[iSlice,j],edge_order=2)
                    
                if validIndxs[-1]+1 - disCtys[-1] > 3:
                    
                    iSlice = np.arange(disCtys[-1],validIndxs[-1]+1)
                    du_dy[iSlice,j] = np.gradient(self.gridVx[iSlice,j],self.gridY[iSlice,j],edge_order=2)
                    
                if len(disCtys) > 1:
                    
                    for disCtyIndx in range(1,len(disCtys)-1):
                        
                        validIndx1 = validIndxs[disCtys[disCtyIndx]+1]
                        validIndx2 = validIndxs[disCtys[disCtyIndx+1]]
                        
                        if validIndx2+1 - validIndx1 > 3:
                            
                            iSlice = np.arange(validIndx1,validIndx2+1)
                            du_dy[iSlice,j] = np.gradient(self.gridVx[iSlice,j],self.gridY[iSlice,j],edge_order=2)
        
        return dv_dx - du_dy
            
    def flipGrid(self,flipAxis=1):
        
        '''flipAxis=1 - horizontal flip
           flipAxis=0 - vertical flip
           flipAxis=None - horizontal and vertical flip'''
        
        destFrame = deepcopy(self)
        
        for gridAttrb in self.frameGridAttrbs[2:]:
            
            setattr(destFrame,gridAttrb,np.flip(getattr(self,gridAttrb),axis=flipAxis))
            
        if flipAxis == 1:
            
            destFrame.gridVx = -destFrame.gridVx
            destFrame.discXc = np.max(destFrame.x) - (np.abs(np.min(destFrame.x)) + destFrame.discXc)
            
        elif flipAxis == 0:
            
            destFrame.gridVy = -destFrame.gridVy
            destFrame.discYc = np.max(destFrame.y) - (np.abs(np.min(destFrame.y)) + destFrame.discYc)
            
        elif flipAxis == None:
            
            destFrame.gridVx = -destFrame.gridVx
            destFrame.gridVy = -destFrame.gridVy
            destFrame.discXc = np.max(destFrame.x) - (np.abs(np.min(destFrame.x)) + destFrame.discXc)
            destFrame.discYc = np.max(destFrame.y) - (np.abs(np.min(destFrame.y)) + destFrame.discYc)
            
        # destFrame.gridWz = destFrame.computeWz2()
        destFrame.flattenData()
        
        return destFrame
    
    def ndimData(self, fstreamRho, fstreamVel):
        
        destFrame = deepcopy(self)
        
        gridAttrbs = deepcopy(self.frameGridAttrbs) #Iterator without Wz
        del gridAttrbs[5]
        
        for gridAttrb in gridAttrbs[:2]:
            
            setattr(destFrame, gridAttrb, getattr(self,gridAttrb)/(self.discDia))
            
        for gridAttrb in gridAttrbs[2:-1]:
            
            setattr(destFrame, gridAttrb, getattr(self,gridAttrb)/fstreamVel)
            
        # destFrame.gridWz = destFrame.computeWz2()
        destFrame.gridWz = destFrame.gridWz*self.discDia*1e-3/fstreamVel#Beware that 1e-3 is used to be compatible with velocity units
        destFrame.flattenData()
        
        destFrame.discXc = destFrame.discXc/(self.discDia)
        destFrame.discYc = destFrame.discYc/(self.discDia)
        destFrame.discDia = 1
        
        return destFrame
    
    def getPolygonIndxs(self,points):
        
        '''Returns grid data containing the truth values at points within defined polygon '''

        polygonObj = Polygon(points)
        
        polygonCheck = np.empty(len(self.x),dtype=bool)
        
        for i,(x, y) in enumerate(zip(self.x, self.y)):
            
            pointObj = Point(x,y)
            polygonCheck[i] = polygonObj.contains(pointObj)
            
        return polygonCheck.reshape(self.gridX.shape)
    
    def getDataLine(self, line, gridAttrb="gridVx"):
        
        return getattr(self,f"{gridAttrb}+Intpr")(line)
        
        
    