#!/usr/bin/env python
"""
Processes a collection of shapefiles to generate distances and coordinates of
points, and then runs a script to calculate an MST using the generated distances
as possible edges.

Alex Zvoleff, aiz2101@columbia.edu
"""

import os
import sys
import getopt

import distGen
import calcTree

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


           # raise Usage("No arguments specified.")
        try:
            #inputFile = argv[1]
            inputFile=r"C:\Users\Selin\Desktop\1000\random1000.shp"
            outputDir =r"C:\Users\Selin\Desktop\1000\random1000"
            algorithm = "Kruskal's Algorithm"
            searchRadius = "100000 meters"
        except IndexError:
            raise Error("Not enough arguments provided to script.")

        distGen.generateDist(inputFile, outputDir, searchRadius)

        calcTree.runAlgorithm(algorithm, outputDir, outputDir)

    except Usage, err:
        print >>sys.stderr, err.msg
        print """
DESCRIPTION:
Given a point shapefile, calculate an MST. This script runs distGen.py and calcTree.py sequentially.

USAGE: fullTreeCalc.py InputShapefile OutputDir Algorithm SearchRadius

SCRIPT PARAMETERS:
InputShapefile: Input shapefile containing a set of points and a \"Weight\" column.
Outputdir: Output directory for network shapefile, nodes.txt, dists.txt and shpProj.prj files.
Algorithm: Either "Kruskal\'s Algorithm" or "Prim\'s Algorithm" (in quotes)
SearchRadius: Radius to limit the size of the resulting set of segments. Units must be provided.

EXAMPLE: fullTreeCalc.py C:\RuhiiraData\Ruhiira.shp C:\RuhiiraData\Output "Kruskal\'s Algorithm" "1000 Meters\""""
        return 2

    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

if __name__ == "__main__":
    sys.exit(main())
