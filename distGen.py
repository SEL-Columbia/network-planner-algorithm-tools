#!/usr/bin/env python
"""
Processes a collection of shapefiles to generate distances between the points,
and coordinates of the points. Search radius is a parameter of the script.

Alex Zvoleff, aiz2101@columbia.edu
"""

# Compute pairwise distances between nodes in a shapefile

import os
import sys
import getopt
import shutil

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

        if len(args) == 0:
            raise Usage("No arguments specified.")

        try:
            inputFile = argv[1]
            outputDir = argv[2]
            searchRadius = argv[3]
        except IndexError:
            raise Error("Not enough arguments provided to script.")

        generateDist(inputFile, outputDir, searchRadius)

    except Usage, err:
        print >>sys.stderr, err.msg
        print """
SCRIPT DESCRIPTION:
Given a point shapefile, generate nodes and distances text files for use in calcTree.py script.

USAGE: distGen.py InputShapefile OutputDir SearchRadius

SCRIPT PARAMETERS:
InputShapefile: Input shapefile containing a set of points and a \"Weight\" column.
Outputdir: Directory in which to output nodes.txt, dists.txt and shpProj.prj datafiles.
SearchRadius: Radius to limit the size of the resulting set of segments. Units must be provided.

EXAMPLE: distGen.py C:\RuhiiraData\Ruhiira.shp C:\RuhiiraData\Output "1000 Meters\""""
        return 2

    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

def generateDist(inputFile, outputPath, searchRadius):
    searchRadius = fileRW.convertToShpUnits(inputFile, searchRadius)
    rootDir, fc = os.path.split(inputFile)
    file, ext = os.path.splitext(fc)

    if not os.path.exists(outputPath):
        try:
            os.mkdir(outputPath)
        except:
            print "ERROR: could not create new directory", outputPath

    print "Loading data from shapefile..."
    nodes = fileRW.readNodesFromShp(inputFile)

    print "Calculating distances..."
    dists = []
   
    nodeCopy = nodes.copy()
    for startNode in nodes.values():
        del nodeCopy[startNode.getID()]
        for endNode in nodeCopy.values():
            dist = ((startNode.getX() - endNode.getX())**2 + (startNode.getY() - endNode.getY())**2)**0.5
         
            if dist < searchRadius:
                dists.append([startNode.getID(), endNode.getID(), dist])

    print "Outputting nodes..."
    nodeFile = outputPath + os.sep + "nodes.txt"
    fileRW.writeNodesToTxt(nodes, nodeFile)

    print "Outputting distances..."
    distFile = outputPath + os.sep + "dists.txt"
    fileRW.writeSegsToTxt(dists, distFile)

    # Copy projection information to shpProj.prj file for later use
    base, ext =  os.path.splitext(fc)
    projFile = os.path.join(rootDir, base + ".prj")
    projDestFile = os.path.join(outputPath, "shpProj.prj")
    try:
        shutil.copy(projFile, projDestFile)
    except:
        print "Could not find projection file"
    return 0

def trimDists(inDistFile, outDirectory, maxDist):
    try:
        inFile = open(inDistFile, "r")
        outFile = open(outDistFile, "r")
    except IOError:
        print "ERROR: could not access distance files"
    for line in inFile:
        node1, node2, weight = line.split()
        if weight <= maxDist:
            outFile.write("%(node1)i %(node2)i %(weight)f\n" %vars())
    inFile.close()
    outFile.close()
    return 0

if __name__ == "__main__":
    sys.exit(main())
