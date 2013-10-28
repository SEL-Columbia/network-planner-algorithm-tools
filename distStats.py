#!/usr/bin/env python
"""
Calculates interhousehold distance statistics for an input shapefile, taking
account of the cumulative population using a weight parameter associated
with each node.

Alex Zvoleff, aiz2101@columbia.edu
"""

import os
import sys
import getopt
import math
import operator

#import pylab
import numpy

from modules import fileRW

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
        for o, a in opts:
            if o in ("-h", "-help"):
                raise Usage('USAGE')
        if len(args) == 0:
            raise Usage("No arguments specified.")

        inputFile = argv[1]

        base, ext =  os.path.splitext(inputFile)
        if ext == '.shp':
            sortField = argv[2] # Only used if inputFile is a shapefile
            costType = argv[3]
            outputPlotfile = argv[4]
            if argv[5] == "#":
                windowLength = None
            else:
                windowLength = int(argv[5])
            if argv[6] == "#":
                title = ""
            else:
                title = argv[6]
            if sortField == "#":
                print "ERROR: sortField must be provided if input file is a shapefile"
                sys.exit(1)
            nodeWeights, segWeights = getStatsFromShapefile(inputFile, sortField)
        elif ext == '.txt':
            costType = argv[2]
            outputPlotfile = argv[3]
            if argv[4] == "#":
                windowLength = None
            else:
                windowLength = int(argv[4])
            if argv[5] == "#":
                title = ""
            else:
                title = argv[5]
            nodeWeights, segWeights = loadStatsFromText(inputFile)
        else:
            print "ERROR: input file must be a shapefile or text file"
            sys.exit(1)

        if costType.lower() == "marginal cost":
            cumNodeWeights=numpy.cumsum(nodeWeights)
            cumSegWeights=numpy.cumsum(segWeights)
            percentConn, costs = calcMargCost(cumNodeWeights, cumSegWeights)
            ylabel = "Marginal distance (m)"
        elif costType.lower() == "mean cost":
            cumNodeWeights=numpy.cumsum(nodeWeights)
            cumSegWeights=numpy.cumsum(segWeights)
            percentConn, costs = calcMeanCost(cumNodeWeights, cumSegWeights)
            ylabel = "Mean distance per connection (m)"
        elif costType.lower() == "total cost":
            cumNodeWeights=numpy.cumsum(nodeWeights)
            cumSegWeights=numpy.cumsum(segWeights)
            percentConn, costs = calcTotalCost(cumNodeWeights, cumSegWeights)
            ylabel = "Total cost (m)"

        if windowLength != None:
            costs = runAvg(costs, windowLength)

        plotStats(percentConn, costs, outputPlotfile, ylabel, title)

    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "For help use --help"
        return 2

    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

def calcMeanCost(cumNodeWeights, cumSegWeights):
    'Returns percent nodes connected versus mean segment weight.'
    percentConn = map(lambda nodeWeight: (nodeWeight / cumNodeWeights[-1]) * 100,
            cumNodeWeights)
    meanSegWeights = map(lambda segWeight, nodeWeight: segWeight / nodeWeight, 
            cumSegWeights, cumNodeWeights)
    return percentConn, meanSegWeights

def calcMeanCostForCluster(cumNodeWeights, cumSegWeights):
    'Returns percent nodes connected versus mean segment weight.'
    numClustersConnected=[]
    for x in range(1,len(cumNodeWeights)+1):
        numClustersConnected.append(x)
    meanSegWeights = map(lambda segWeight,nodeWeight: (float(segWeight) / nodeWeight), 
            cumSegWeights, cumNodeWeights)

 
    percentConn = map(lambda transformerWeight: (float(transformerWeight) / numClustersConnected[-1]) * 100,
            numClustersConnected)
    
    return numClustersConnected,percentConn, meanSegWeights

def calcMargCost(cumNodeWeights, cumSegWeights):
    'Returns percent nodes connected versus marginal segment weight.'
    percentConn = map(lambda nodeWeight: (nodeWeight / cumNodeWeights[-1]) * 100,
            cumNodeWeights)
    # Remove 100 percent as there is no meaningful marg cost for 100 percent:
    percentConn = percentConn[:-1] 
    margSegWeights = map(lambda nextWeight, thisWeight: nextWeight - thisWeight,
            cumSegWeights[1:], cumSegWeights[:-1])
    # Now divide margSegWeights so that units are margSegWeight / nodeWeight
    margSegWeights = map(lambda nextNodeWeight, thisNodeWeight, margSegWeight:
            margSegWeight / (nextNodeWeight - thisNodeWeight),
            cumNodeWeights[1:], cumNodeWeights[:-1], margSegWeights)
    return percentConn, margSegWeights

def calcTotalCost(cumNodeWeights, cumSegWeights):
    'Returns percent nodes connected versus total segment weight.'
    percentConn = map(lambda nodeWeight: (nodeWeight / cumNodeWeights[-1]) * 100,
            cumNodeWeights)
    return percentConn, cumSegWeights

def calcNearestNeighborIndex(shapefile, area=None):
    nodes = fileRW.readNodesFromShp(shapefile)
    # Calculate avg nearest neighbor distance:
    minDists = {}
    nodeCopy = nodes.copy()
    minXY = [None, None]
    maxXY = [None, None]
    for startNode in nodes.values():
        startX = startNode.getX()
        startY = startNode.getY()
        # Store minimum coordinates to make a bounding rectangle
        if startX < minXY[0] or minXY[0] == None:
            minXY[0] = startX
        elif startX > maxXY[0] or maxXY[0] == None:
            maxXY[0] = startX
        if startY < minXY[1] or minXY[1] == None:
            minXY[1] = startY
        elif startY > maxXY[1] or maxXY[1]  == None:
            maxXY[1] = startY
        startNodeID = startNode.getID()
        del nodeCopy[startNode.getID()]
        for endNode in nodeCopy.values():
            endNodeID = endNode.getID()
            dist = ((startX - endNode.getX())**2 + 
                    (startY - endNode.getY())**2)**(.5)
            if not minDists.has_key(endNodeID) or dist < minDists[endNodeID]:
                minDists[endNodeID] = dist
            if not minDists.has_key(startNodeID) or dist < minDists[startNodeID]:
                minDists[startNodeID] = dist
    if area == None:
        area = (maxXY[0] - minXY[0]) * (maxXY[1] - minXY[1])
    numPts = len(minDists.values())
    avgNeighDist = sum(minDists.values())/numPts
    expNeighDist = .5 * math.sqrt(area / numPts)
    print 'Average distance to nearest neighbor:', avgNeighDist
    print 'Expected distance to nearest neighbor:', expNeighDist
    print 'Nearest neighbor index:', avgNeighDist / expNeighDist
    return avgNeighDist / expNeighDist

def runAvg(series, windowLength):
    'Calculate a moving average using a rectangular window to smooth cost plot.'
    window = numpy.ones(windowLength, 'd')
    runAvgSeries = numpy.convolve(window/window.sum(), series, mode = 'same')
    # Address endpoint problem by using an asymmetric window.
    for n in xrange(1,windowLength/2):
        runAvgSeries[-n] = numpy.mean(series[-(n+windowLength/2):])
        runAvgSeries[n-1] = numpy.mean(series[:windowLength/2+n])
    return runAvgSeries

def getStatsFromShapefile(shapefile, sortField):
    'Gets segment weight and population information from a network shapefile.'
    net = fileRW.readNetFromShp(shapefile)
    segDict = net.getEdgesDict()

    # Determine order of values based on sortField
    if sortField == "FID":
        sortedSegs = segDict.values()
    else:
        sortFieldValues = fileRW.readFieldFromShp(shapefile, sortField)
        sortedItems = sortFieldValues.items()
        sortedItems.sort(key=operator.itemgetter(1))
        sortedSegs = []
        for item in sortedItems:
            FID = item[0]
            sortedSegs.append(segDict[FID])
    segWeights = []
    nodeWeights = []
    nodesInNet = set([])
    for seg in sortedSegs:
        node1, node2 = seg.getNode1(), seg.getNode2()
        endPts = [node1, node2]
        newNodeWeights = [node1.getWeight(), node2.getWeight()]
        length = seg.getWeight()
        newNodeWeightSum = float(0)
        for endPt, newNodeWeight in zip(endPts, newNodeWeights):
            if endPt not in nodesInNet:
                nodesInNet.add(endPt)
                newNodeWeightSum += newNodeWeight
        nodeWeights.append(newNodeWeightSum)
        segWeights.append(length)
    return nodeWeights, segWeights

def getStatsFromTree(tree, sortByWeight=False):
    segs = tree.getEdges()
    if sortByWeight == True:
        segs.sort()
    segWeights = []
    nodeWeights = []
    nodesInNet = set([])
    for seg in segs:
        nodes = seg.getNodes()
        newNodeWeightSum = float(0)
        for node in nodes:
            if node not in nodesInNet:
                nodesInNet.add(node)
                newNodeWeightSum += node.getWeight()
        nodeWeights.append(newNodeWeightSum)
        segWeights.append(seg.getWeight())
    return nodeWeights, segWeights

def loadStatsFromText(statsFile):
    'Loads statistics from a text file.'
    try:
        inFile = open(statsFile, "r")
    except:
        print "ERROR: Could not load statistics file."
        return 1
    x = []
    y = []
    for line in inFile:
        xi, yi = line.split()
        x.append(float(xi))
        y.append(float(yi))
    inFile.close()
    return x, y

def loadStatsFromTextAsDict(statsFile):
    'Loads statistics from a text file.'
    try:
        inFile = open(statsFile, "r")
    except:
        print "ERROR: Could not load Dict file."
        return 1
    Dict = {}
    for line in inFile:
        xi, yi = line.split()
        Dict[float(xi)]=float(yi)
    inFile.close()
    return Dict
    



def writeStatsToText(statsFile, x, y):
    'Writes statistics to a text file for later reuse.'
    outFile = open(statsFile,"w")
    for xi, yi in zip(x, y):
            outFile.write("%(xi)f %(yi)f\n" %vars())
    outFile.close()
    return 0

def loadKeyFromText(keyFile):
    'Loads statistics key from a text file.'
    try:
        inFile = open(keyFile, "r")
    except:
        print "ERROR: Could not load statistics file."
        return 1
    costsKey = {}
    for line in inFile:
        percent, startFID = line.split()
        costsKey[percent] = startFID
    inFile.close()
    return costsKey

def writeKeyToText(keyFile, costsKey):
    'Writes statistics key file to text for later reuse.'
    outFile = open(keyFile,"w")
    for percent in sorted(costsKey.keys()):
        startFID = costsKey[percent]
        outFile.write("%(percent)f %(startFID)i\n" %vars())
    outFile.close()
    return 0

def plotStats(x, y, outputFile, ylabel, title=""):
    'Uses pylab to produce a plot of cumulative statistics.'
    pylab.plot(x, y, 'k-', linewidth=2)
    pylab.xlabel("Households connected (%)", fontsize='large')
    pylab.ylabel(ylabel, fontsize='large')
    pylab.xlim(0, 100)
    pylab.xticks([0, 20, 40, 60, 80, 100], ('0%', '20%', '40%', '60%', '80%',
        '100%'), fontsize='large')
    pylab.yticks(fontsize='large')
    pylab.title(title, fontsize='x-large')
    pylab.savefig(outputFile, dpi=300)
    pylab.clf()
    return 0

def regress(x, y, order):
    n = len(y)
    coeff = numpy.polyfit(x, y, order)
    yr = numpy.polyval(coeff, x)
    error = math.sqrt(sum((yr-y)**2)/n)
    return (coeff, error, yr)

def plotSegHist(shapeFile):
    'Plot histogram of network segment lengths.'
    segs, nodes = fileRW.getNetFromShape
    #TODO: Finish coding function.

def thinStats(stats, length):
    'Used to thin a list. Useful for printing small plots.'
    if len(stats) < length:
        print "ERROR: original length less than shortened length"
    dx = len(stats) / length
    shortStats = []
    for n in xrange(len(stats)):
        if n % dx == 0:
            shortStats.append(stats[n])
    return shortStats

def roundStats(percentConn, costs, statsLength):
    """Rounds percents as necessary so as to fit percents and costs to a
    lower resolution (so that there are 100/n values)."""
    dx = 100. / statsLength
    costsDict = {}
    indices = {}
    counter = 0
    for percent, cost in zip(percentConn, costs):
        remainder = percent % dx
        if remainder > .5 * dx:
            percent = percent - remainder + dx
        else:
            percent = percent - remainder
        percent = numpy.round(percent, 3) # Round to avoid diffs due to machine precision
        if costsDict.has_key(percent) and costsDict[percent] < cost:
            pass
        else:
            costsDict[percent] = cost
            indices[percent] = counter
        counter += 1
    percentConn = []
    costs = []
    for percent in sorted(costsDict.keys()):
        percentConn.append(percent)
        costs.append(costsDict[percent])
    return percentConn, costs, indices

if __name__ == "__main__":
    sys.exit(main())
