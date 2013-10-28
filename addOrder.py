#!/usr/bin/env python
"""
Uses prim's algorithm to build a network successively from a given point.
Codes the segments according to the order they are added, and according
to their percentile in the network.

Alex Zvoleff, aiz2101@columbia.edu
"""

import os
import sys
import getopt

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

        inputShapefile = argv[1] # Precalculated MST
        startPt = int(argv[2])
        fieldName = argv[3]
        maxNodes = argv[4]

        print "Reading network from input shapefile..."
        tree = fileRW.readNetFromShp(inputShapefile)
        nodes = tree.getNodes()
        segs = tree.getEdges()

        print "Calculating tree..."
        if maxNodes == "#":
            numNodes = len(nodes)
        else:
            numNodes = float(maxNodes)
        tree = prims.primsAlg(segs, numNodes, startPt)

        print "Writing order to shapefile..."
        # segOrder is a dictionary where keys are seg FIDs,
        # and values are the order added to the network
        segOrder = {}
        n = 0
        for seg in tree.getEdges():
           segOrder[seg.getID()] = n
           n += 1
        fileRW.writeFieldToShp(inputShapefile, segOrder, fieldName)
        
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "For help use --help"
        return 2

    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

if __name__ == "__main__":
    sys.exit(main())
