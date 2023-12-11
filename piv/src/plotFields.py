from matplotlib import pyplot as plt
import numpy as np
from numpy.ma import masked_array
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors
import matplotlib.patches as patches

from src.ioUtilsNoStitch import caseData
from src.miscTools import axesBboxSize

class plotFields:
    
    def __init__(self, caseDataPath, axBounds=None, figObjExt=None): #axBounds syntax [left, right, top, bottom]
            
        phases = np.arange(18, 360, 36)
        
        self.leftAx = None
        self.rightAx = None
        self.topAx = None
        self.bottomAx = None
        
        if axBounds != None:
            
            self.leftAx = axBounds[0]
            self.rightAx = axBounds[1]
            self.topAx = axBounds[2]
            self.bottomAx = axBounds[3]
        
        self.data = caseData(caseDataPath)
        self.vMax = self.data.fstreamVelc + 0.5
        
        if self.data.phaseQty == 1:
            
            if figObjExt != None:
                
                self.figVel = figObjExt[0]
                self.figVort = figObjExt[1]
                
                if self.data.caseNo == 6:
                    
                    self.axVel = self.figVel.axes[0]
                    self.axVort = self.figVort.axes[0]
                    
                elif self.data.caseNo == 7:
                    
                    self.axVel = self.figVel.axes[1]
                    self.axVort = self.figVort.axes[1]
                    
            else:
            
                self.figVel, self.axVel = plt.subplots(layout="tight")
                self.figVort, self.axVort = plt.subplots(layout="tight")
            
            figData = getattr(self.data,f"phase0").combinedFrames
            velPlot, vortPlot = self.genFieldPlots(figData, self.axVel, self.axVort)
            
            if figObjExt == None:
                
                divider = make_axes_locatable(self.axVel)
                caxObj = divider.append_axes("right", size="2%", pad=0.05)
                cbar = self.figVel.colorbar(velPlot, cax=caxObj)
                cbar.set_label(r"$\|V\|/V_\infty\ [-]$")
                
                self.axVel.set_xlabel("x/D [-]")
                self.axVel.set_ylabel("y/D [-]")
                axesBboxSize(0.01 * 1.02 * (self.rightAx - self.leftAx), 0.01 * (self.topAx - self.bottomAx), self.axVel)
                
                divider = make_axes_locatable(self.axVort)
                caxObj = divider.append_axes("right", size="2%", pad=0.05)
                cbar = self.figVort.colorbar(vortPlot, cax=caxObj)
                cbar.set_label(r"$\omega_z\ \times D/V_\infty [-]$")
                
                self.axVort.set_xlabel("x/D [-]")
                self.axVort.set_ylabel("y/D [-]")
                axesBboxSize(0.01 * 1.02 * (self.rightAx - self.leftAx), 0.01 * (self.topAx - self.bottomAx), self.axVort)
            
        else:
            
            self.figVel, self.axVel = plt.subplots(5, 2, figsize=(2 * 0.01 * 1.02 * (self.rightAx - self.leftAx), 5 * 0.01 * (self.topAx - self.bottomAx)), sharex=True, sharey=True, layout="constrained")
            self.figVort, self.axVort = plt.subplots(5, 2, figsize=(2 * 0.01 * 1.02 * (self.rightAx - self.leftAx), 5 * 0.01 * (self.topAx - self.bottomAx)), sharex=True, sharey=True, layout="constrained")
            
            phaseIndex = 0
            
            for i in range(5):
                
                for j in range(2):
                    
                    figData = getattr(self.data,f"phase{phaseIndex}").combinedFrames
                    self.genFieldPlots(figData, self.axVel[i][j], self.axVort[i][j])
                    self.axVel[i][j].set_title(f"phase = {phases[phaseIndex]}"+r"$^\degree$")
                    self.axVort[i][j].set_title(f"phase = {phases[phaseIndex]}"+r"$^\degree$")
                    
                    phaseIndex = phaseIndex + 1
                    
            cbarVel = self.figVel.colorbar(self.velPlot, ax=self.axVel[1:4, 1], aspect=100, pad=0)
            cbarVel.set_label(r"$\|V\|/V_\infty\ [-]$")
            
            cbarVort = self.figVort.colorbar(self.vortPlot, ax=self.axVort[1:4, 1], aspect=100, pad=0)
            cbarVort.set_label(r"$\omega_z\ \times D/V_\infty [-]$")
            
            for axVel in self.axVel[:,0]:
                
                axVel.set_ylabel("y/D [-]")
                
            for axVort in self.axVort[:,0]:
                
                axVort.set_ylabel("y/D [-]")
                
            
            for axVel in self.axVel[-1,:]:
                
                axVel.set_xlabel("x/D [-]")
                
            for axVort in self.axVort[-1,:]:
                
                axVort.set_xlabel("x/D [-]")
                
    def genFieldPlots(self, figData, axVelObj, axVortObj):
        
        maskedX = masked_array(figData.gridPosX, mask=~figData.gridIsValid)
        maskedY = masked_array(figData.gridPosY, mask=~figData.gridIsValid)
        maskedVelX = masked_array(figData.gridVelX, mask=~figData.gridIsValid)
        maskedVelY = masked_array(figData.gridVelY, mask=~figData.gridIsValid)
        maskedVel = masked_array(figData.gridVel, mask=~figData.gridIsValid)/figData.phaseVinf
        maskedVort = masked_array(figData.gridVortZ, mask=~figData.gridIsValid)*(self.data.discDia*10**-3/figData.phaseVinf)
        
        print("Min vel: ", np.min(maskedVel))
        print("Max vel: ", np.max(maskedVel))
        
        cdict = {'red': [(0.0, 33/255, 33/255),
                         (1/6, 103/255, 103/255),
                         (2/6, 209/255, 209/255),
                         (3/6, 247/255, 247/255),
                         (4/6, 253/255, 253/255),
                         (5/6, 239/255, 239/255),
                         (1.0, 178/255, 178/255)],
                 'green': [(0.0, 102/255, 102/255),
                           (1/6, 169/255, 169/255),
                           (2/6, 229/255, 229/255),
                           (3/6, 247/255, 247/255),
                           (4/6, 219/255, 219/255),
                           (5/6, 138/255, 138/255),
                           (1.0, 24/255, 24/255)],
                 'blue': [(0.0, 172/255, 172/255),
                          (1/6, 207/255, 207/255),
                          (2/6, 240/255, 240/255),
                          (3/6, 247/255, 247/255),
                          (4/6, 199/255, 199/255),
                          (5/6, 98/255, 98/255),
                          (1.0, 43/255, 43/255)]}
        
        cdict2 = {'red' : [(0.0, 0.0, 0.0),
                           (0.49999, 0.0, 1.0),
                           (0.5, 1.0, 100/255),
                           (0.50001, 100/255, 100/255),
                           (1.0, 1.0, 1.0)],
                  'green' : [(0.0, 0.0, 0.0),
                             (0.49999, 0.0, 1.0),
                             (0.5, 1.0, 0.0),
                             (0.50001, 0.0, 0.0),
                             (1.0, 0.0, 0.0)],
                  'blue' : [(0.0, 1.0, 1.0),
                            (0.49999, 100/255, 1.0),
                            (0.5, 1.0, 100/255),
                            (0.50001, 100/255, 0.0),
                            (1.0, 0.0, 0.0)],}
        
        davidCmap = colors.LinearSegmentedColormap("davidCmap", segmentdata=cdict)
        kiranCmap = colors.LinearSegmentedColormap("kiranCmap", segmentdata=cdict2)
        
        divnorm = colors.TwoSlopeNorm(vmin=np.min(maskedVel), vcenter=1.0, vmax=np.max(maskedVel))
        
        velPlot = axVelObj.contourf(figData.gridPosX/self.data.discDia, figData.gridPosY/self.data.discDia, maskedVel, levels=np.linspace(0.5, np.max(maskedVel),25), cmap="seismic", extend="min", norm=divnorm)
        axVelObj.quiver(maskedX[::24, ::24]/self.data.discDia, maskedY[::24, ::24]/self.data.discDia, maskedVelX[::24, ::24], maskedVelY[::24, ::24])
        
        vortPlot = axVortObj.contourf(figData.gridPosX/self.data.discDia, figData.gridPosY/self.data.discDia, maskedVort, levels=np.linspace(-75*0.0835,75*0.0835,50), cmap="seismic", extend="both")
        axVortObj.quiver(maskedX[::24, ::24]/self.data.discDia, maskedY[::24, ::24]/self.data.discDia, maskedVelX[::24, ::24], maskedVelY[::24, ::24])
        
        self.velPlot = velPlot
        self.vortPlot = vortPlot
        
        if self.leftAx != None:
            
            axVelObj.set_xlim(left=self.leftAx/self.data.discDia, right=self.rightAx/self.data.discDia)
            axVelObj.set_ylim(bottom =self.bottomAx/self.data.discDia, top=self.topAx/self.data.discDia)
            
            axVortObj.set_xlim(left=self.leftAx/self.data.discDia, right=self.rightAx/self.data.discDia)
            axVortObj.set_ylim(bottom =self.bottomAx/self.data.discDia, top=self.topAx/self.data.discDia)
            
        xmin, xmax = axVelObj.get_xlim()
        ymin, ymax = axVelObj.get_ylim()
        xy = (xmin,ymin)
        width = xmax - xmin
        height = ymax - ymin
        
        #Patches
        hatchedPatchVel = patches.Rectangle(xy, width, height, fc = "dimgrey", zorder=-10)
        hatchedPatchVort = patches.Rectangle(xy, width, height, fc = "dimgrey", zorder=-10)
        axVelObj.add_patch(hatchedPatchVel)
        axVortObj.add_patch(hatchedPatchVort)
            
        return velPlot, vortPlot
    
    def save(self, dirPath):
        
        self.figVel.savefig(dirPath+f"{self.data.case}"+"VelField.png", dpi=300, pad_inches=0)
        self.figVort.savefig(dirPath+f"{self.data.case}"+"VortField.png", dpi=300, pad_inches=0)
            
            