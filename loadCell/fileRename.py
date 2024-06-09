import os
import glob
import sys
from pathlib import Path

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

fileList = sorted(glob.glob(dataDir + "/loadCellTest/20NME/*RPM[1-9][0-9][0-9].txt"))

for filePath in fileList:
    
    pathObj = Path(filePath)
    newFilePath = pathObj.stem[:-3] + "0" + pathObj.stem[-3:] + pathObj.suffix
    # pathObj.rename(Path(pathObj.parent,newFilePath))
    
    
    
    
    
