
#!/usr/bin/env python
"""
Ayse Selin Kocaman
"""

import os
import sys
import getopt
import shutil
import time
import pickle

import pylab
import numpy

import distStats

from modules import prims
from modules import fileRW
from modules import network

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Error(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.GetoptError, msg:
            raise Usage(msg)

        #if len(args) == 0:
            #raise Usage("No arguments specified.")

        startTime = time.time()    

        # argv[1] must be an input shapefile containing a network shapefile 
        # created with the calcTree.py script. This shapefile must contains a
        # precalculated spanning tree with columns containing the node IDs of
        # the endpoints of each segment ("pt1" and "pt2") and the weights of
        # these endpoints, "pt1Weight" and "pt2weight".
        try:
            #inputShapefile = argv[1] # a pre-calculated MST
            inputShapefile=r"C:\Users\Selin\Desktop\real\tiby_all\TibyQB2005Pts\MV.shp"
            #outputPath = argv[2] # Path for the output of the batchPrims script
            outputPath=r"C:\Users\Selin\Desktop\real\tiby_all\TibyQB2005Pts"
            #locationTitle = argv[3] # Name of the site, for example "Sauri, Kenya"
            locationTitle="Tiby"
        
        except IndexError:
            raise Error("Not enough arguments provided to script.")


        MVtree=fileRW.readNetFromShp(inputShapefile)
        batchPrims(MVtree,outputPath,locationTitle)


        
    except Usage, err:
        print >>sys.stderr, err.msg
        return 2
    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

    

def batchPrims(tree,outputPath,locationTitle):
    dataPath = outputPath + os.sep + "data" 
    bestPath = outputPath + os.sep + "best"
    if not os.path.exists(outputPath):
        os.mkdir(outputPath)
    if not os.path.exists(dataPath):
        os.mkdir(dataPath)
    if not os.path.exists(bestPath):
        os.mkdir(bestPath)
    startFIDs = []
    resolution = 500 #Alex needed this.
    segs = tree.getEdges()
    nodes = tree.getNodes()
    nodesDict=tree.getNodeDict()

    nodeSegmentsDict = prims.buildAssocDict(segs)
    for node in nodes:
        startFIDs.append(node.getID())
    n = 1
    numStartPts = len(startFIDs)
    maxNodes = len(nodes)
    for startFID in startFIDs:
            if n % 10 == 0:
                print "Run %i of %i..." %(n, numStartPts)

            tree = prims.primsAlg(segs, maxNodes, startFID, [], nodeSegmentsDict)
            #fileRW.genShapefile(tree, outputPath + ".prj", outputPath + os.sep + "prim%i.shp" %(startFID))
            statsFile = dataPath + os.sep + "batchPrims%i.txt" %(startFID)
            nodeWeights, totalCost,MV,LV = getStatsFromTree(tree,nodesDict,outputPath,startFID)# bunun yerine yeni bisi yaz....
            cumNodeWeights = numpy.cumsum(nodeWeights) # amount of population/demand that is connected to the grid
            cumTotalCost = numpy.cumsum(totalCost) # cost or sum of benefit function
            cumMV=numpy.cumsum(MV)
            cumLV=numpy.cumsum(LV)
            writeStatsToTextforClusterBased(statsFile,cumNodeWeights,cumTotalCost,cumMV, cumLV)
            n += 1
    print "***Done with Prims calculations***"
    #print "Elapsed time: %s" %(elapsedTime(startTime))

    bestCostsByPercent,bestCostsByNumber,MVByNumber,MVByPercent,LVByNumber,LVByPercent,numClustersConn,percentConn = calcBestStatsClusterBased(dataPath)
   
    statsFile = bestPath + os.sep + "batchPrimsBestByNumber.txt" # specifies lowest mean cost/highest benefit for each penetration rate
    writeStatsToTextforClusterBased(statsFile, numClustersConn, bestCostsByNumber,MVByNumber,LVByNumber)
    statsFile = bestPath + os.sep + "batchPrimsBestByPercent.txt"
    writeStatsToTextforClusterBased(statsFile, percentConn, bestCostsByPercent,MVByPercent,LVByPercent)

    
    statsPlotByNumber = bestPath + os.sep + "batchPrimsBestByNumber.pdf"
    statsPlotByPercent = bestPath + os.sep + "batchPrimsBestByPercent.pdf"
    #plotTitle="Mbola"
    plotStatsForClusterBased(numClustersConn, bestCostsByNumber,MVByNumber,LVByNumber, statsPlotByNumber, "Mean Cost", locationTitle)
    plotStatsForClusterBased(percentConn, bestCostsByPercent,MVByPercent,LVByPercent, statsPlotByPercent, "Mean Cost", locationTitle)



def getStatsFromTree(tree,centers,outputDir,startFID,sortByWeight=False):
    segs = tree.getEdges()
    DictFile= outputDir + os.sep + "LVCostDict.txt"
    LVCostDict=distStats.loadStatsFromTextAsDict(DictFile)
    
    if sortByWeight == True:
        segs.sort()
    totalCost = []
    MVLineMeter=[]
    LVLineMeter=[]
    nodeWeights = []
    
    totalCost.append(LVCostDict[startFID]+5000)
    LVLineMeter.append(float(LVCostDict[startFID])/10)
    MVLineMeter.append(0)
    nodeWeights.append(centers[startFID].getWeight())
    
    nodesInNet = set([])
    nodesInNet.add(centers[startFID])
    for seg in segs:
        nodes = seg.getNodes()
        newNodeWeightSum = float(0)
        for node in nodes:
            if node not in nodesInNet:
                nodesInNet.add(node)
                newNodeWeightSum += node.getWeight()
                newLV=LVCostDict[node.getID()]
        nodeWeights.append(newNodeWeightSum)
        totalCost.append(seg.getWeight()*25+newLV+5000) #5000 additional transformer cost
        MVLineMeter.append(seg.getWeight())
        LVLineMeter.append(float(newLV)/10)
    return nodeWeights,totalCost,MVLineMeter,LVLineMeter

def calcBestStatsClusterBased(dataPath):
    """Calculates best costs for each percentage out of costs calculated from a
    large set of starting points. Resolution is the resolution the data is
    sampled over, in percent."""
    bestCostsByPercent = {}
    bestCostsByNumber={}
    MVByNumber={}#MV Line Length in meter
    LVByNumber={} #LV Line Length in meter
    MVByPercent={}
    LVByPercent={}
    bestCostsArrayByPercent=[]# Stores the best costs for each percentage penetration
    bestCostsArrayByNumber=[]
    MVByNumberArray=[]
    LVByNumberArray=[]
    MVByPercentArray=[]
    LVByPercentArray=[]
    
    for file in os.listdir(dataPath):
        if not file.startswith('batchPrims') or not file.endswith('.txt'):
            continue
        startFID = int(file[10:-4]) # Remove 'batchPrims' prefix and '.txt' suffix
        filePath = os.path.join(dataPath,file)
        cumNodeWeights,cumSegWeights,cumMVLength,cumLVLength = loadStatsFromTextforClusterBased(filePath)
       
        numClustersConnected,percentConnClusters, costs,meanMV,meanLV = calcMeanCostForCluster(cumNodeWeights, cumSegWeights,cumMVLength,cumLVLength)
        #print "percentConn", percentConn
        for number,cost,MV,LV in zip(numClustersConnected, costs, meanMV, meanLV):
            #print "number, cost", number, cost
            if bestCostsByNumber.has_key(number) and (bestCostsByNumber[number] < cost):
                pass
            else:
                bestCostsByNumber[number] = cost
                LVByNumber[number]=LV
                MVByNumber[number]=MV
                
        for percent, cost ,MV,LV in zip(percentConnClusters, costs, meanMV, meanLV):
            #print "percent, cost", percent, cost
            if bestCostsByPercent.has_key(percent) and (bestCostsByPercent[percent] < cost):
                pass
            else:
                bestCostsByPercent[percent] = cost
                LVByPercent[percent]=LV
                MVByPercent[percent]=MV

    for percent in sorted(bestCostsByPercent.keys()):
         bestCostsArrayByPercent.append(bestCostsByPercent[percent])
    for percent in sorted(LVByPercent.keys()):
         LVByPercentArray.append(LVByPercent[percent])
    for percent in sorted(MVByPercent.keys()):
         MVByPercentArray.append(MVByPercent[percent])

         
    for number in sorted(bestCostsByNumber.keys()):
         bestCostsArrayByNumber.append(bestCostsByNumber[number])
    for number in sorted(LVByNumber.keys()):
         LVByNumberArray.append(LVByNumber[number])
    for number in sorted(MVByNumber.keys()):
         MVByNumberArray.append(MVByNumber[number]) 
                
    #print  bestCostsByPercent,bestCostsByNumber   
    return bestCostsArrayByPercent,bestCostsArrayByNumber,MVByNumberArray,MVByPercentArray,LVByNumberArray,LVByPercentArray,numClustersConnected,percentConnClusters

def calcMeanCostForCluster(cumNodeWeights, cumSegWeights,cumMVLength,cumLVLength):
    'Returns percent nodes connected versus mean segment weight.'
    numClustersConnected=[]
    for x in range(1,len(cumNodeWeights)+1):
        numClustersConnected.append(x)
    meanSegWeights = map(lambda segWeight,nodeWeight: (float(segWeight) / nodeWeight), 
            cumSegWeights, cumNodeWeights)
    
    meanMV = map(lambda MVLength,nodeWeight: (float(MVLength) / nodeWeight), 
            cumMVLength, cumNodeWeights)
    
    meanLV = map(lambda LVLength,nodeWeight: (float(LVLength) / nodeWeight), 
            cumLVLength, cumNodeWeights)
    
    percentConnClusters = map(lambda transformerWeight: (float(transformerWeight) / numClustersConnected[-1]) * 100,
            numClustersConnected)

    return numClustersConnected,percentConnClusters, meanSegWeights,meanMV,meanLV

def writeStatsToTextforClusterBased(statsFile,x,y,z,t):
    'Writes statistics to a text file for later reuse.'
    outFile = open(statsFile,"w")
    for xi, yi ,zi,ti in zip(x, y,z,t):
            outFile.write("%(xi)f %(yi)f %(zi)f %(ti)f \n" %vars())
    outFile.close()
    return 0

def loadStatsFromTextforClusterBased(statsFile):
    'Loads statistics from a text file.'
    try:
        inFile = open(statsFile, "r")
    except:
        print "ERROR: Could not load statistics file."
        return 1
    x = []
    y = []
    z = []
    t = []
    for line in inFile:
        xi, yi,zi,ti = line.split()
        x.append(float(xi))
        y.append(float(yi))
        z.append(float(zi))
        t.append(float(ti))
    inFile.close()
    return x, y, z, t



def plotStatsForClusterBased(x, y1,y2,y3 ,outputFile, ylabel, title=""):
    'Uses pylab to produce a plot of cumulative statistics.'
    xlen=len(y1)
    
    pylab.ylabel('Mean LV and Mv Lengths(meter)', fontsize='x-large')
    pylab.xlabel("Number of transformers connected", fontsize='x-large')
    MV=pylab.plot(x, y2, 'b--', linewidth=5)
    LV=pylab.plot(x, y3, 'y.-', linewidth=5)
    pylab.ylim(0,100)
    pylab.twinx()
    meanCost=pylab.plot(x, y1, 'm-', linewidth=5)
    if x[0]==1:
        pylab.xlabel("Number of transformers connected", fontsize='x-large')
        ticks=[]
        for i in range(0,xlen/20+1):
            ticks.append(i*20)
        ticks.append(xlen)
        if title=="Ruhiira":
            pylab.xticks(ticks, ('0/%i'%(xlen), '20/%i'%(xlen), '40/%i'%(xlen),'49/%i'%(xlen) ), fontsize='xlarge')
        
        if title=="Pampaida":
            pylab.xticks(ticks, ('0/%i'%(xlen), '20/%i'%(xlen), '40/%i'%(xlen), '60/%i'%(xlen), '80/%i'%(xlen),'100/%i'%(xlen),'115/%i'%(xlen)), fontsize='x-large') 
        if title=="Mbola":
            pylab.xticks(ticks, ('0/%i'%(xlen), '20/%i'%(xlen), '40/%i'%(xlen), '60/%i'%(xlen), '80/%i'%(xlen),'100/%i'%(xlen),'120/%i' %(xlen),'139/%i' %(xlen)), fontsize='large')
        if title=="Tiby": #tick leri de degistir
            pylab.xticks(ticks, ('0/%i'%(xlen), '5/%i'%(xlen), '10/%i'%(xlen),'15/%i'%(xlen)),fontsize='large')

    else:
        pylab.xlabel("Transformers connected (%)", fontsize='xlarge')
        pylab.xticks([0, 20, 40, 60, 80, 100], ('0%', '20%', '40%', '60%', '80%',
        '100%'), fontsize='x-large')
        pylab.xlim(0,100)
        
    pylab.ylabel(ylabel, fontsize='x-large')
    pylab.ylim(min(min(y1),min(y2),min(y3)), max(max(y1),max(y2),max(y3))+100)
    pylab.xlabel("Number of transformers connected", fontsize='x-large')
    pylab.yticks(fontsize='large')
    pylab.title(title, fontsize='x-large')
    pylab.legend((MV,LV,meanCost),('MeanMV','MeanLV','MeanCost'),loc="upper left")
    pylab.savefig(outputFile, dpi=300)
    pylab.clf()
    return 0

def makePastaPlot(dataPath, plotFile, plotTitle):
    pylab.figure()
    for file in os.listdir(dataPath):
        if ('batchPrims' not in file) and ('.txt' not in file):
            continue
        startFID = int(file[10:][:-4]) #remove 'batchPrims' prefix and '.txt' suffix
        filePath = os.path.join(dataPath,file)
        nodeWeights, segWeights = distStats.loadStatsFromText(filePath)
        percentConn, costs = distStats.calcMeanCost(nodeWeights, segWeights)
        pylab.plot(percentConn, costs, 'k-', linewidth=.5)
    pylab.xlabel("Households connected (%)", fontsize='large')
    pylab.ylabel("Mean distance per connection (m)", fontsize='large')
    pylab.xlim(0, 100)
    pylab.xticks([0, 20, 40, 60, 80, 100], ('0%', '20%', '40%', '60%', '80%',
        '100%'), fontsize='large')
    pylab.yticks(fontsize='large')
    pylab.title(plotTitle, fontsize='x-large')
    if os.path.exists(plotFile):
        os.remove(plotFile)
    pylab.savefig(plotFile, dpi=300)

def elapsedTime(startTime):
    elapsed = int(time.time() - startTime)
    hours = elapsed / 3600
    minutes = (elapsed - hours * 3600) / 60
    seconds = elapsed - hours * 3600 - minutes * 60
    return "%ih %im %is" %(hours, minutes, seconds)

if __name__ == "__main__":
    sys.exit(main())
