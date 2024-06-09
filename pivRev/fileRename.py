import os
import sys
import glob
from pathlib import Path

#Configuring relative file locations
homeDir = os.path.dirname(__file__)
srcDir = os.path.join(homeDir,"src")
plotDir = os.path.join(homeDir,"plots")
dataDir = os.path.join(homeDir,"data")

# fileList = glob.glob(dataDir + "/*Case[0-9]*.dat")

# Renaming tecData mean files
# for filePath in fileList:
    
#     pathObj = Path(filePath)
    
#     if len(pathObj.stem.split("_")[1]) == 1:
#         newFilePath = pathObj.stem[:3] + "Mean" + pathObj.stem[3:7] + "0" + pathObj.stem[7] + "Phase" + "0" + pathObj.stem.split("_")[1] + ".dat"
        
#     else:
        
#         newFilePath = pathObj.stem[:3] + "Mean" + pathObj.stem[3:7] + "0" + pathObj.stem[7] + "Phase"  + pathObj.stem.split("_")[1] + ".dat"
        
    # pathObj.rename(Path(pathObj.parent,newFilePath))
    
#Renaming log files
fileList = glob.glob(dataDir + "/*Log*Phase[0-9].txt")

for filePath in fileList:
    
    pathObj = Path(filePath)
    
    # newFilePath = "p70LogCase05" + pathObj.stem.capitalize() + ".txt"
    
    newFilePath = pathObj.stem[:-1] + "0" + pathObj.stem[-1] + ".txt"
    
    # pathObj.rename(Path(pathObj.parent,newFilePath))
    
    