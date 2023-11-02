import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from miscTools import animFrameUpdate

class cyclePlot:
    
    def __init__(self, windObj, noWindObj, aeroObj, staticObj):
        
        self.rawAnimFig, self.rawAnimAxs = plt.subplots()
        self.filtAnimFig, self.filtAnimAxs = plt.subplots()
        
        self.windObj = windObj
        self.noWindObj = noWindObj
        self.aeroObj = aeroObj
        self.staticObj = staticObj
        
    def genAnim(self, plotDir):
        
        #Initalize raw animation
        self.rawAnimAxs.plot(self.aeroObj.cycleTime[0], self.windObj.cycleForce[0], color='b', label="total")
        self.rawAnimAxs.plot(self.aeroObj.cycleTime[0], self.noWindObj.cycleSyncForce[0], color='r', label="inertial")
        self.rawAnimAxs.plot(self.aeroObj.cycleTime[0], self.aeroObj.cycleForce[0], color='g', label="aero")
        self.rawAnimAxs.plot(self.aeroObj.cycleTime[0], self.windObj.cycleAvgForce, color='b', linestyle="--", label="total")
        self.rawAnimAxs.plot(self.aeroObj.cycleTime[0], self.noWindObj.cycleAvgSyncForce, color='r', linestyle="--", label="inertial")
        self.rawAnimAxs.plot(self.aeroObj.cycleTime[0], self.aeroObj.cycleAvgForce, color='g', linestyle="--", label="aero")
        self.rawAnimAxs.set_xlabel("Time [s]")
        self.rawAnimAxs.set_ylabel(f"$C_T$ [-]")
        self.rawAnimAxs.legend()
        
        #Initialize filt animation
        self.filtAnimAxs.plot(self.aeroObj.cycleTime[0], self.windObj.cycleFiltForce[0], color='b', label="total")
        self.filtAnimAxs.plot(self.aeroObj.cycleTime[0], self.noWindObj.cycleSyncFiltForce[0], color='r', label="inertial")
        self.filtAnimAxs.plot(self.aeroObj.cycleTime[0], self.aeroObj.cycleFiltForce[0], color='g', label="aero")
        self.filtAnimAxs.plot(self.aeroObj.cycleTime[0], self.windObj.cycleAvgFiltForce, color='b', linestyle="--", label="total")
        self.filtAnimAxs.plot(self.aeroObj.cycleTime[0], self.noWindObj.cycleAvgSyncFiltForce, color='r', linestyle="--", label="inertial")
        self.filtAnimAxs.plot(self.aeroObj.cycleTime[0], self.aeroObj.cycleAvgFiltForce, color='g', linestyle="--", label="aero")
        self.filtAnimAxs.set_xlabel("Time [s]")
        self.filtAnimAxs.set_ylabel(f"$C_T$ [-]")
        self.filtAnimAxs.legend()
        
        #Save animations
        rawAnim = animation.FuncAnimation(fig=self.rawAnimFig, func=animFrameUpdate, frames=self.aeroObj.trimSurgeCycles, fargs=(self,False))
        filtAnim = animation.FuncAnimation(fig=self.filtAnimFig, func=animFrameUpdate, frames=self.aeroObj.trimSurgeCycles, fargs=(self,True))
        
        rawAnim.save(filename=plotDir+"/rawAnim"+self.windObj.uid+".gif", writer="pillow")
        filtAnim.save(filename=plotDir+"/filtAnim"+self.windObj.uid+".gif", writer="pillow")
        
        return rawAnim, filtAnim
        
    def genPlot(self, cycle, plotDir):
        
        rawFig, rawAxs = plt.subplots()
        filtFig,filtAxs = plt.subplots()
        
        xMin = self.aeroObj.cycleNdimTime[0]
        xMax = self.aeroObj.cycleNdimTime[-1]        
        
        #rawAxs.plot(self.aeroObj.cycleNdimTime, self.windObj.cycleForce[cycle], color='b', label="total")
        #rawAxs.plot(self.aeroObj.cycleNdimTime, self.noWindObj.cycleSyncForce[cycle], color='r', label="inertial")
        #rawAxs.plot(self.aeroObj.cycleNdimTime, self.aeroObj.cycleForce[cycle], color='g', label="aero")
        #rawAxs.plot(self.aeroObj.cycleNdimTime, self.windObj.cycleAvgForce, color='b', linestyle="--", label="total")
        #rawAxs.plot(self.aeroObj.cycleNdimTime, self.noWindObj.cycleAvgSyncForce, color='r', linestyle="--", label="inertial")
        rawAxs.plot(self.aeroObj.cycleNdimTime, self.aeroObj.cycleAvgForce, color='g', linestyle="--", label="aero")
        rawAxs.fill_between(self.aeroObj.cycleNdimTime, self.aeroObj.cycleAvgForce - self.aeroObj.cycleStdForce, self.aeroObj.cycleAvgForce + self.aeroObj.cycleStdForce, alpha=0.2)
        rawAxs.set_xlabel("Time [-]")
        rawAxs.set_ylabel(f"$C_T$ [-]")
        rawAxs.legend()
        
        #filtAxs.plot(self.aeroObj.cycleNdimTime, self.windObj.cycleFiltForce[cycle], color='b', label="total")
        #filtAxs.plot(self.aeroObj.cycleNdimTime, self.noWindObj.cycleSyncFiltForce[cycle], color='r', label="inertial")
        #filtAxs.plot(self.aeroObj.cycleNdimTime, self.aeroObj.cycleFiltForce[cycle], color='g', label="aero")
        #filtAxs.plot(self.aeroObj.cycleNdimTime, self.windObj.cycleAvgFiltForce, color='b', linestyle="--", label="total")
        #filtAxs.plot(self.aeroObj.cycleNdimTime, self.noWindObj.cycleAvgSyncFiltForce, color='r', linestyle="--", label="inertial")
        #filtAxs.plot(self.aeroObj.cycleNdimTime, self.aeroObj.cycleAvgForce, color='r', label="surgeAeroRaw")    
        filtAxs.plot(self.aeroObj.cycleNdimTime, self.aeroObj.cycleAvgFiltForce, color='g', label="surgeAero")
        #filtAxs.fill_between(self.aeroObj.cycleNdimTime, self.aeroObj.cycleAvgFiltForce - self.aeroObj.cycleStdForce, self.aeroObj.cycleAvgFiltForce + self.aeroObj.cycleStdForce, alpha=0.2)
        filtAxs.hlines(np.mean(self.staticObj.forceSeries), xmin = xMin, xmax = xMax, color='k', label="staticAero")
        filtAxs.set_xlim(xMin, xMax)
        filtAxs.set_xlabel("Time [-]")
        filtAxs.set_ylabel(f"$C_T$ [-]")
        filtAxs.legend()
        filtAxs.grid()
        
        #rawFig.savefig(fname=plotDir+f"/rawPlotCycle{cycle}"+self.windObj.uid+".png")
        filtFig.savefig(fname=plotDir+f"/filtAvgThrust"+self.windObj.uid+".png")
        
        
        