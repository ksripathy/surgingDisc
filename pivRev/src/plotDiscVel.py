from matplotlib import pyplot as plt
import numpy as np
from matplotlib import colors
import pickle

class plotDiscVel:

    def __init__(self, dataDir, plotDir):
        
        self.dataDir = dataDir
        self.plotDir = plotDir
    
    def loadCase(self, discPor, caseID):
        
        with open(self.dataDir + f"/p{discPor}Case{caseID}.pickle","rb") as handle:
            dynCaseObj = pickle.load(handle)
            
        if int(caseID) < 3:
        
            with open(self.dataDir + f"/p{discPor}Case06.pickle","rb") as handle:
                staticCaseObj = pickle.load(handle)
                
        else:
            
             with open(self.dataDir + f"/p{discPor}Case07.pickle","rb") as handle:
                staticCaseObj = pickle.load(handle)
            
        return dynCaseObj, staticCaseObj
    
    def staticPlot(self, discPor, discVelAttrb="discVxcExtp"):
        
        staticCase6 = self.loadCase(discPor,"06")[0]
        staticCase7 = self.loadCase(discPor,"07")[0]
        
        plotX = np.linspace(-0.5,0.5,21)
        plotY1 = getattr(staticCase6,discVelAttrb)
        # plotY2 = getattr(staticCase6,discVelAttrb)
        plotY3 = getattr(staticCase7,discVelAttrb)
        # plotY4 = getattr(staticCase7,discVelAttrb)
        
        staticFig, staticAx = plt.subplots()
        staticAx.plot(plotX,plotY1,marker="x",label="Case6")
        # staticAx.plot(plotX,plotY2,marker="x",label="Case6,Linear")
        staticAx.plot(plotX,plotY3,marker="x",label="Case7")
        # staticAx.plot(plotX,plotY4,marker="x",label="Case8,RBF",color="tab:orange",linestyle="--")
        
        staticAx.set_xlabel("y/D [-]")
        staticAx.set_ylabel(r"$V_D/V_{\infty} [-]$")
        
        staticAx.legend()
        
        staticFig.savefig(self.plotDir + f"/p{discPor}{discVelAttrb}Static.png",dpi=600,pad_inches=0,bbox_inches="tight")
        
    
    def pcolorPlot(self, discPor, caseID, discVelAttrb="discVxcExtp"):
        
        dynCaseObj, staticCaseObj = self.loadCase(discPor,caseID)
        
        staticDiscVel = getattr(staticCaseObj,discVelAttrb)[10]
            
        plotPhases = np.linspace(-0.1*np.pi,2.1*np.pi,12)
        axisPhases = np.linspace(0,2*np.pi,11)
        plotSpanLocs = np.linspace(-0.55,0.55,12)
        axisSpanLocs = np.linspace(-0.5,0.5,11)
        
        plotX, plotY = np.meshgrid(plotPhases,plotSpanLocs)
        plotZ = np.zeros((11,11))
        plotZ[:,:-1] = np.roll(getattr(dynCaseObj,discVelAttrb)[::2,:10],-2,axis=1) - staticDiscVel #roll for init phase at 90degree phase
        plotZ[:,-1] = plotZ[:,0]
        
        divnorm = colors.TwoSlopeNorm(vmin=np.min(plotZ), vcenter=0, vmax=np.max(plotZ))
        
        pcolorFig, pcolorAx = plt.subplots()
        pcolorArtist = pcolorAx.pcolor(plotX, plotY, plotZ, edgecolor="k",linewidths=1, cmap="seismic", norm=divnorm)
        pcolorAx.set_xticks(axisPhases)
        pcolorAx.set_yticks(axisSpanLocs)
        pcolorAx.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 5))
        pcolorAx.xaxis.set_major_formatter(Multiple(5).formatter())
        pcolorAx.set_xlabel(r"$\phi [rad]$")
        pcolorAx.set_ylabel(r"$y/D [-]$")
        pcolorCbar = pcolorFig.colorbar(pcolorArtist)
        pcolorCbar.set_label(r"$(V_D - V_{D,static})/V_{\infty} [-]$")
        
        pcolorFig.savefig(self.plotDir + f"/p{discPor}Case{caseID}{discVelAttrb}Pcolor.png",dpi=600,pad_inches=0,bbox_inches="tight")
        
    def casesPlot(self, discPor, caseIDArr = ["00", "01", "02", "03", "04", "05"], spanID="10", discVelAttrb="discVxcExtp"):
        
        dynCaseObjArr = []
        staticCaseObjArr = []
        staticDiscVelArr = np.zeros(len(caseIDArr))
        
        lineColors = ["tab:blue", "tab:blue", "tab:orange", "tab:green", "tab:green", "tab:red"]
        lineStyles = ["-", "--", "-", "-", "--", "-"]
        lineMarkers = ["x", "^", "x", "x", "^", "x"]
        
        for caseID in caseIDArr:
            
            dynCaseObj, staticCaseObj = self.loadCase(discPor,caseID)
            dynCaseObjArr.append(dynCaseObj)
            staticCaseObjArr.append(staticCaseObj)
            
        plotX = np.linspace(0*np.pi,2*np.pi,11)
        axisPhases = np.linspace(0,2*np.pi,11)
        
        casesFig, casesAx = plt.subplots()
        spanNo = int(spanID)
        
        for i,(dynCaseObj, staticCaseObj) in enumerate(zip(dynCaseObjArr, staticCaseObjArr)):
            
            staticDiscVelArr[i] = getattr(staticCaseObj,discVelAttrb)[spanNo]
            
            plotY = np.zeros(11)
            discVelArr = np.roll(getattr(dynCaseObj,discVelAttrb)[spanNo,:10],-2)
            plotY[:-1] = discVelArr
            plotY[-1] = discVelArr[0]
            
            casesAx.plot(plotX,plotY, color=lineColors[i], marker=lineMarkers[i], markerfacecolor = "none", linestyle=lineStyles[i],label=dynCaseObj.uid[3:])
        
        staticDiscVel = np.mean(staticDiscVelArr)
        
        casesAx.axhline(staticDiscVel, label="Static",linestyle = "--", color="k")
        casesAx.axvline(0.5*np.pi,color="green",linewidth=1,linestyle="dotted")
        casesAx.axvline(1.5*np.pi,color="red",linewidth=1,linestyle="dotted")
            
        casesAx.set_xticks(axisPhases)
        casesAx.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 5))
        casesAx.xaxis.set_major_formatter(Multiple(5).formatter())
        casesAx.set_xlabel(r"$\phi [rad]$")
        casesAx.set_ylabel(r"$V_D/V_{\infty} [-]$")
        casesAx.legend()
        
        casesFig.savefig(self.plotDir + f"/p{discPor}{discVelAttrb}Span{spanID}Cases.png",dpi=600,pad_inches=0,bbox_inches="tight")
        
    def dynInflowPlot(self, discPor, caseIDArr = ["00", "01", "02", "03", "04", "05"], spanID="10", discVelAttrb="discVxcExtp"):
        
        dynCaseObjArr = []
        staticCaseObjArr = []
        staticDiscVelArr = np.zeros(len(caseIDArr))
        
        lineColors = ["tab:blue", "tab:blue", "tab:orange", "tab:green", "tab:green", "tab:red"]
        lineStyles = ["-", "--", "-", "-", "--", "-"]
        lineMarkers = ["x", "^", "x", "x", "^", "x"]
        
        for caseID in caseIDArr:
            
            dynCaseObj, staticCaseObj = self.loadCase(discPor,caseID)
            dynCaseObjArr.append(dynCaseObj)
            staticCaseObjArr.append(staticCaseObj)
            
        phases = np.linspace(0*np.pi,2*np.pi,11)
        axisPhases = np.linspace(0,2*np.pi,11)
        plotX = np.cos(phases + 0.5*np.pi) #Additional 90deg since the plots are made with it as init phase. Theoretical surge motion velocity
            
        dynFig, dynAx = plt.subplots()
        spanNo = int(spanID)
        
        for i,(dynCaseObj, staticCaseObj) in enumerate(zip(dynCaseObjArr, staticCaseObjArr)):
            
            staticDiscVelArr[i] = getattr(staticCaseObj,discVelAttrb)[spanNo]
            
            plotY = np.zeros(11)
            discVelArr = np.roll(getattr(dynCaseObj,discVelAttrb)[spanNo,:10],-2)
            plotY[:-1] = discVelArr
            plotY[-1] = discVelArr[0]
            
            dynAx.plot(plotX,plotY,color=lineColors[i], marker=lineMarkers[i], markerfacecolor = "none", linestyle=lineStyles[i],label=dynCaseObj.uid[3:])
            
        staticDiscVel = np.mean(staticDiscVelArr)
        
        dynAx.plot(0,staticDiscVel)
        dynAx.axvline(-1,color="green",linewidth=1,linestyle="dotted")
        dynAx.axvline(1,color="red",linewidth=1,linestyle="dotted")
        
        dynAx.set_xticks(plotX)
        dynAx.set_xlabel(r"$u/V_{\infty} [-]$")
        dynAx.set_ylabel(r"$V_D/V_{\infty} [-]$")
        dynAx.legend()
        
        dynFig.savefig(self.plotDir + f"/p{discPor}{discVelAttrb}Span{spanID}dynInflow.png",dpi=600,pad_inches=0,bbox_inches="tight")
        
    def fieldPlot(self, discPor, caseID):
        
        dynCaseObj, staticCaseObj = self.loadCase(discPor,caseID)
        
        phaseObjArr = []
        
        for i,phaseName in enumerate(dynCaseObj.phaseObjNames):
            
            phaseObj = getattr(dynCaseObj,phaseName)
            phaseObjArr.append(phaseObj)
            
            velFig, velAx = plt.subplots()
            vortFig, vortAx = plt.subplots()
            
            X = dynCaseObj.gridX - phaseObj.discXc
            Y = dynCaseObj.gridY - phaseObj.discYc
            V = phaseObj.gridV
            Vx = phaseObj.gridVx
            Vy = phaseObj.gridVy
            Wz = phaseObj.gridWz
            
            divnorm = colors.TwoSlopeNorm(vmin=np.min(V), vcenter=1.0, vmax=np.max(V))
            mask = np.ma.getmaskarray(V)
            maskData = np.zeros(X.shape)
            maskData[mask] = 1
            discMask = np.ma.masked_array(maskData,np.logical_not(mask))
            
            trim = 50
            
            velPlot = velAx.contourf(X[:,trim:-trim],Y[:,trim:-trim],Vx[:,trim:-trim],cmap="seismic",norm=divnorm, levels=np.linspace(np.min(Vx),1.5,21),extend="max")
            velAx.contourf(X[:,trim:-trim],Y[:,trim:-trim],discMask[:,trim:-trim],cmap="binary")
            velAx.quiver(X[::24,trim:-trim:24],Y[::24,trim:-trim:24],Vx[::24,trim:-trim:24],Vy[::24,trim:-trim:24])
            velAx.set_xlabel("x/D [-]")
            velAx.set_ylabel("y/D [-]")
            velCbar = velFig.colorbar(velPlot, orientation="horizontal")
            velCbar.set_label(r"$Vx/V_{\infty} [-]$")
            velAx.set_aspect("equal")
            
            velFig.savefig(self.plotDir + f"/p{discPor}Case{caseID}{phaseName}VelField.png",dpi=600,pad_inches=0,bbox_inches="tight")
            
            
            vortPlot = vortAx.contourf(X[:,trim:-trim],Y[:,trim:-trim],Wz[:,trim:-trim],cmap="seismic",levels=np.linspace(-100,100,21),extend="both")
            vortAx.contourf(X[:,trim:-trim],Y[:,trim:-trim],discMask[:,trim:-trim],cmap="binary")
            vortAx.quiver(X[::24,trim:-trim:24],Y[::24,trim:-trim:24],Vx[::24,trim:-trim:24],Vy[::24,trim:-trim:24])
            vortAx.set_xlabel("x/D [-]")
            vortAx.set_ylabel("y/D [-]")
            vortCbar = velFig.colorbar(vortPlot, orientation="horizontal")
            vortCbar.set_label(r"$\omega_z D/V_\infty [-]$")
            vortAx.set_aspect("equal")
            
            vortFig.savefig(self.plotDir + f"/p{discPor}Case{caseID}{phaseName}VortField.png",dpi=600,pad_inches=0,bbox_inches="tight")
            
            
            
            
            
            
        
        
        
        
        
        
        
        
#Class for plotting axis ticks as multiple of pi
#https://stackoverflow.com/questions/40642061/how-to-set-axis-ticks-in-multiples-of-pi-python-matplotlib
class Multiple:
    def __init__(self, denominator=2, number=np.pi, latex='\pi'):
        self.denominator = denominator
        self.number = number
        self.latex = latex
        
    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)
    
    def formatter(self):
        return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))
    
def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
    
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    
    def _multiple_formatter(x, pos):
        den = denominator
        num = int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            if num==1:
                return r'$%s$'%latex
            elif num==-1:
                return r'$-%s$'%latex
            else:
                return r'$%s%s$'%(num,latex)
        else:
            if num==1:
                return r'$\frac{%s}{%s}$'%(latex,den)
            elif num==-1:
                return r'$\frac{-%s}{%s}$'%(latex,den)
            else:
                return r'$\frac{%s%s}{%s}$'%(num,latex,den)
            
    return _multiple_formatter
        