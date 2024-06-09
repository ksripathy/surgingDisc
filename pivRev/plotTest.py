import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

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

class Multiple:
    def __init__(self, denominator=2, number=np.pi, latex='\pi'):
        self.denominator = denominator
        self.number = number
        self.latex = latex
        
    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)
    
    def formatter(self):
        return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

with open(dataDir + "/pickledData/p45Case02.pickle","rb") as handle:
    caseObj = pickle.load(handle)

pcolorX = np.linspace(0,2*np.pi,11)
pcolorY = np.linspace(-0.55,0.55,12)

surgeCumPhases = np.linspace(0.1*np.pi,2.3*np.pi,12)
surgePhases = np.array(surgeCumPhases)
surgePhases[10:] = surgePhases[10:] - (2*np.pi)
radialLocs = np.linspace(-0.5,0.5,21)

gridPcolorY,gridPcolorX = np.meshgrid(pcolorX,pcolorY)

divnorm = colors.TwoSlopeNorm(vmin=np.min(caseObj.indVelRBF.Vxc[::2,:10]-0.57), vcenter=0, vmax=np.max(caseObj.indVelRBF.Vxc[::2,:10]-0.57))

fig,ax = plt.subplots()
pcolorPlot = ax.pcolor(gridPcolorY,gridPcolorX,np.roll(caseObj.indVelRBF.Vxc[::2,:10]-0.57,-2,axis=1),edgecolor="k",linewidths=1, cmap="seismic", norm=divnorm)
plt.xticks(pcolorX)
plt.yticks(radialLocs[::2])
# ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 5))
# ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
ax.xaxis.set_minor_formatter(Multiple(5).formatter())
ax.set_xlabel(r"$\phi [rad]$")
ax.set_ylabel(r"$y/D [-]$")
pcolorBar = fig.colorbar(pcolorPlot)
pcolorBar.set_label(r"$V_D/V_{\infty} [-]$")