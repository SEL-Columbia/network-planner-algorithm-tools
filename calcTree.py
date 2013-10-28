#!/usr/bin/env python
"""
Runs networking algorithm on a precalculated set of segments and node
coordinates. Can use either Kruskal's or Prim's algorithms.

Alex Zvoleff, aiz2102@columbia.edu
"""

import os
import sys
import getopt

from modules import fileRW
from modules import kruskals
from modules import prims

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

        try:
            inputDir = argv[1]
            outputDir = argv[2]
            algorithm = argv[3]
            if argv[4] == "#":
                cutoff = None
            else:
                cutoff = float(argv[4])
        except IndexError:
            raise Error("Not enough arguments provided to script.")

        runAlgorithm(algorithm, inputDir, outputDir, cutoff)

    except Usage, err:
        print >>sys.stderr, err.msg
        print """
SCRIPT DESCRIPTION:
Calculate a network given input nodes, distance, and projection text files.

USAGE: calcTree.py InputDir OutputDir Algorithm Cutoff

SCRIPT PARAMETERS:
InputDir: Input directory containing output of distGen.py script
OutputDir: Directory in which to output network.shp shapefile
Algorithm: Either "Kruskal\'s Algorithm" or "Prim\'s Algorithm" (in quotes)
Cutoff: Optional parameter to further limit the search radius. To leave blank use \"#\"

EXAMPLE: calcTree.py C:\RuhiiraData C:\RuhiiraData\Network "Kruskal\'s Algorithm" "1000\""""
        return 2

    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

def runAlgorithm(algorithm, inputDir, outputDir, cutoff=None):
    if not os.path.exists(outputDir):
        try:
            os.mkdir(outputDir)
        except:
            raise Error("Could not create new directory " + inputFile)
    distFile = inputDir + os.sep + "dists.txt"
    nodeFile = inputDir + os.sep + "nodes.txt"
    projFile = inputDir + os.sep + "shpProj.prj"
    shapePath = outputDir + os.sep + "network.shp"
    if not os.path.exists(inputDir):
        print inputDir + " does not exist"
        return 1
    elif not os.path.exists(outputDir): 
        print outputDir + " does not exist"
        return 1
    elif not os.path.exists(distFile): 
        print distFile + " does not exist"
        return 1
    elif not os.path.exists(nodeFile): 
        print nodeFile + " does not exist"
        return 1
    elif not os.path.exists(projFile): 
        print projFile + " does not exist"
        return 1

    print "Generating network hehehh:"

    nodes = fileRW.readNodesFromTxt(nodeFile)
    testSegments = fileRW.readSegsFromTxt(distFile, nodes, cutoff)

    numNodes = len(nodes)
    print numNodes
   
    
    print "\tCalculating tree..."
    if algorithm == "Kruskal's Algorithm":
        tree = kruskals.kruskalsAlg(testSegments, numNodes)
    elif algorithm == "Prim's Algorithm":
        tree = prims.primsAlg(testSegments, numNodes,firstNodeID=0)
    else:
        print "Specified algorithm does not exist"

    print "\tNodes:\t\t" + str(tree.numNodes())
    print "\tConnections:\t" + str(tree.numEdges())
    print "\tSubnet(s):\t\t" + str(tree.numSubnets())

    # Alex says sum benefit function over tree to get total benefit for tree
    # Can compare trees using this to pick best from set of trees
    print "\tTotal Weight:\t" + str(tree.getTotalEdgeWeight())

    if tree.numSubnets() == 1:
        MID = tree.getTotalEdgeWeight() / tree.getTotalNodeWeight()
        print "\tMID (100):\t\t" + str(MID)
    else:
        print "\tMultiple subnets. No MID calculated."

    print "Generating network shapefile..."
    fileRWReturnCode = fileRW.genShapefile(tree, projFile, shapePath)
    if fileRWReturnCode == 0:
        print "\tNetwork saved as shapefile at: " + shapePath
        return 0
    else:
       return 1

if __name__ == "__main__":
    sys.exit(main())
