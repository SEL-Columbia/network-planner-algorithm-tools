#!/usr/bin/env python
"""
Implements the "Composite Prim's" algorithm to estimate how optimal network cost
varies with penetration rate. Given a precalculated spanning network as a
shapefile, the script tests an input set of starting points and compiles the
results, returning the optimal cost for each penetration rate, and the node ID
of the optimal starting point for each penetration rate.


Alex Zvoleff, aiz2101@columbia.edu
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

        if len(args) == 0:
            raise Usage("No arguments specified.")

        startTime = time.time()    

        # argv[1] must be an input shapefile containing a network shapefile 
        # created with the calcTree.py script. This shapefile must contains a
        # precalculated spanning tree with columns containing the node IDs of
        # the endpoints of each segment ("pt1" and "pt2") and the weights of
        # these endpoints, "pt1Weight" and "pt2weight".
        try:
            inputShapefile = argv[1] # a pre-calculated MST

            outputPath = argv[2] # Path for the output of the batchPrims script
            locationTitle = argv[3] # Name of the site, for example "Sauri, Kenya"

            if argv[4] == "#":
                startPtFile = None
            else:
                startPtFile = argv[4]
            if argv[5] == "#":
                # sum of node weights
                maxNodes = None
            else:
                maxNodes = int(argv[5])
            if argv[6] == "#":
                searchRadius = None
            else:
                searchRadius = argv[6]
        except IndexError:
            raise Error("Not enough arguments provided to script.")


        dataPath = outputPath + os.sep + "data" 
        bestPath = outputPath + os.sep + "best"
        if not os.path.exists(outputPath):
            os.mkdir(outputPath)
        if not os.path.exists(dataPath):
            os.mkdir(dataPath)
        if not os.path.exists(bestPath):
            os.mkdir(bestPath)

        length = 500 # Resolution used in calcBestStats function

        # This code reads in the precalculated spanning network from
        # inputShapefile to create a set of nodes and a set of segments
        # connecting these nodes ("edges") for use in Prim's algorithm.
        print "Reading distances & nodes..."
        tree = fileRW.readNetFromShp(inputShapefile)
        segs = tree.getEdges()
        #print segs
        nodes = tree.getNodes()
        if maxNodes == None:
            maxNodes = len(nodes)

        print "Building node dictionary..."
        #for efficiency, build nodeDict so it can be passed to Prim's and reused
        nodeDict = prims.buildAssocDict(segs)

        # Network starting points to be tested by the algorithm can be input in
        # one of three ways:
        #     (1) If startPtFile is not specified, then all the nodes in the
        #     network will be tested as potential starting points.
        #
        #     (2) If startPtFile is a point shapefile with "Weight" and "FID"
        #     columns, then a list of starting nodes will be created by making a
        #     list of, for each node in startPtFile, the closest node in
        #     inputShapefile.
        #    
        #     (3) If startPtFile is a text file containing one node ID per
        #     line, the corresponding nodes in inputShapefile will be used as
        #     network starting points.


        # Note that the weight here is not using when building the network
        # It is used at the end and might be used for resolving time-different images
        print "Setting up starting points..."
        startFIDs = []
        if startPtFile != None:
            base, ext =  os.path.splitext(startPtFile)
            if ext == '.shp':
                # Make more error tolerant if shapefile has weights or not
                startNodes = fileRW.readNodesFromShp(startPtFile)
                for startNode in startNodes:
                    # so nearestNode is used in case the specified start nodes are not actually exactly nodes in the tree
                    # so we get the closest one for each startnode.
                    startFIDs.append(nearestNode(startNode, nodes))
            elif ext == '.txt':
                inFile = open(startPtFile, 'r')
                for line in inFile:
                    startFIDs.append(int(line))
        else:
            for node in nodes:
                startFIDs.append(node.getID())

        # The code for the actual batch Prim's ("Composite Prim's") begins
        # below. Networks are built from each starting point using Prim's
        # algorithm, and the results saved for later analysis.
        n = 1
        numStartPts = len(startFIDs)
        print "Elapsed time: %s" %(elapsedTime(startTime))
        print "***Running Prim's Algorithm***"
        print "Run %i of %i..." %(n, numStartPts)
        
        for startFID in startFIDs:
            if n % 10 == 0:
                print "Run %i of %i..." %(n, numStartPts)

            # so for electricity maxNodes was always everything go for full tree
            # but for healthcare they used maxNodes to limit the support of the tree (limit the number of structures or say limit supply of electricity)
            tree = prims.primsAlg(segs, maxNodes, startFID, [], nodeDict)

            # statsFile stores the total length (or cumulative weight) of the
            # edges, and the total weight (or population connected to the grid)
            # for each penetration rate for a network built from each
            # particular starting point.
            
            statsFile = dataPath + os.sep + "batchPrims%i.txt" %(startFID)
            nodeWeights, segWeights = distStats.getStatsFromTree(tree)

            #print "nodeweight", nodeWeights
            cumNodeWeights = numpy.cumsum(nodeWeights) # amount of population/demand that is connected to the grid
            #print "cumNode", cumNodeWeights
            cumSegWeights = numpy.cumsum(segWeights) # cost or sum of benefit function
            distStats.writeStatsToText(statsFile, cumNodeWeights, cumSegWeights)

            n += 1
            if n % 10 == 0:
                print "Elapsed time: %s" %(elapsedTime(startTime))

        print "***Done with Prims calculations***"
        print "Elapsed time: %s" %(elapsedTime(startTime))

        # This code processes the statsFiles produced from each network starting
        # point to built a composite picture of the "optimal" network and
        # starting point for each penetration rate. "statsKeyFile" stores the
        # nodeID of the best starting point for each penetration rate. Here
        # "statsFile" is a text file storing, in rows, the total node weight and total
        # segment weight of the optimal network for each penetration rate.
        print "Finding optimal costs..."
        # Use length to bin (and scale cumulatively) the cumulative demand into penetration rate 
        # and cumulative benefit into mean cost
        bestNodeWeights, bestSegWeights, bestCostsKey = calcBestStats(dataPath, length)
        statsFile = bestPath + os.sep + "batchPrimsBest.txt" # specifies lowest mean cost/highest benefit for each penetration rate
        distStats.writeStatsToText(statsFile, bestNodeWeights, bestSegWeights)
        statsKeyFile = bestPath + os.sep + "batchPrimsBestKey.txt" # specifies best starting point for each penetration rate 
        distStats.writeKeyToText(statsKeyFile, bestCostsKey)
        percentConn, meanCosts = distStats.calcMeanCost(bestNodeWeights,
                bestSegWeights)
        statsPlot = bestPath + os.sep + "batchPrimsBest.pdf"
        plotTitle = "%s Composite Prim's" %(locationTitle)
        distStats.plotStats(percentConn, meanCosts, statsPlot, "MID (m)", plotTitle)
        print "Elapsed time: %s" %(elapsedTime(startTime))

        # The "spaghetti plot" plots the network weight versus penetration rate
        # for networks starting from all of the input starting points. Creating
        # this plot for large networks or a large number of starting points
        # takes a lot of time and a lot of ram, and may crash the script, so
        # this is the final action the script takes, and it is commented out by
        # default.
        #print "Making spaghetti plot..."
        #pastaPlotFile = bestPath + os.sep + "batchPrimsPasta.png"
        #pastaPlotTitle = "%s Spaghetti Plot" %(locationTitle)
        #makePastaPlot(dataPath, pastaPlotFile, pastaPlotTitle)

        print "Total time: %s" %(elapsedTime(startTime))
        
    except Usage, err:
        print >>sys.stderr, err.msg
        print """
SCRIPT DESCRIPTION:
Implements the "Composite Prim's" algorithm to estimate how optimal network
cost varies with penetration rate. Given a precalculated spanning network as a
shapefile, the script tests an input set of starting points and compiles the
results, returning the optimal cost for each penetration rate, and the node ID
of the starting point necessary to acheive that cost.

USAGE: batchPrims.py NetworkShapefile OutputDir Location StartPtFile
    \ MaxNodes SearchRadius

SCRIPT PARAMETERS:
NetworkShapefile: Input shapefile containing a precalculated minumim spanning tree.
Outputdir: Directory in which to output results of the script.
Location: Title (e.g. name of site) to be used in plots.
StartPtFile: See note below.
MaxNodes: Maxmimum number of nodes to include in each network.
SearchRadius: Radius to limit the size of the resulting set of segments. Units must be provided.

NOTE ON StartPtFile:
Network starting points to be tested by the algorithm can be input in one of
three ways:
    (1) If startPtFile is not specified, (entered as "#") then all the nodes in
    the network will be tested as potential starting points.

    (2) If startPtFile is a point shapefile with "Weight" and "FID" columns,
    then a list of starting nodes will be created by making a list of, for each
    node in startPtFile, the closest node in inputShapefile.
   
    (3) If startPtFile is a text file containing one node ID per line, the
    corresponding nodes in inputShapefile will be used as network starting
    points.

EXAMPLES:
(1) To use all nodes as starting points with no search radius or max nodes:
batchPrims.py C:\RuhiiraData\\network.shp C:\RuhiiraData\\batchPrims 
    \ "Ruhiira" "#" "#" "#"\"

(2) To limit the number of nodes to 1000 and use a search radius of 1000 meters,
while using a shapefile of potential starting points:
batchPrims.py C:\RuhiiraData\\network.shp C:\RuhiiraData\\batchPrims 
    \ "Ruhiira" "C:\RuhiiraData\\startPts.shp" "1000" "1000 Meters"\""""
        return 2

    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

def nearestNode(startNode, nodes):
    'Calculates the nearest node in the \'nodes\' set to a given startNode.'
    minDist = None
    startPtNum = 0
    for node in nodes.values():
        dist = ((startNode.getX() - node.getX())**2 + 
                (startNode.getY() - node.getY())**2)**(.5)
        if (dist < minDist) | (minDist == None):
            minDist = dist
            closestNode = node.getID()
    return closestNode

def calcBestStats(dataPath, resolution):
    """Calculates best costs for each percentage out of costs calculated from a
    large set of starting points. Resolution is the resolution the data is
    sampled over, in percent."""
    bestCosts = {} # Stores the best costs for each percentage penetration
    bestCostsKey = {} # Stores the FIDs of the best starting points
    weights = {}
    for file in os.listdir(dataPath):
        if not file.startswith('batchPrims') or not file.endswith('.txt'):
            continue
        startFID = int(file[10:-4]) # Remove 'batchPrims' prefix and '.txt' suffix
        filePath = os.path.join(dataPath,file)
        cumNodeWeights, cumSegWeights = distStats.loadStatsFromText(filePath)
        percentConn, costs = distStats.calcMeanCost(cumNodeWeights, cumSegWeights)
        percentConn, costs, indices = distStats.roundStats(percentConn, costs, resolution)
        for percent, cost in zip(percentConn, costs):
            if bestCosts.has_key(percent) and (bestCosts[percent] < cost):
                pass
            else:
                bestCosts[percent] = cost
                bestCostsKey[percent] = startFID
                # Find the cumulative node/seg weights assoc. with this percent
                weightIndex = indices[percent]
                weights[percent] = (cumNodeWeights[weightIndex], cumSegWeights[weightIndex])
    bestNodeWeights = []
    bestSegWeights = []
    for percent in sorted(weights.keys()):
        nodeWeight, segWeight = weights[percent]
        bestNodeWeights.append(nodeWeight)
        bestSegWeights.append(segWeight)
    return bestNodeWeights, bestSegWeights, bestCostsKey

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
