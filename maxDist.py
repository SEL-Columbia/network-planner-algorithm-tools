#!/usr/bin/env python
"""
Given an input point shapefile, calculates the maximum distance between any
point and its nearest neighbor, and returns that distance.

Alex Zvoleff, aiz2101@columbia.edu
"""
import sys
import os
import getopt

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
            raise Usage("No arguments specified")

        inputFile = sys.argv[1]

        calcMaxDist(inputFile)

    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "For help use --help"
        return 2

    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

def calcMaxDist(inputFile):
    """Calculates the maximum distance of any node in the inputfile to its nearest
    neighbor"""
    nodes = fileRW.readNodesFromShp(inputFile)
    coords = []
    for node in nodes.values():
        coords.append([node.getX(), node.getY()])

    maxDist = None
    startPtNum = 0
    for startPt in coords:
        thisNearest = None
        # Loop over the coords with the startPt cut out.
        endPtIndices = range(0,startPtNum) + \
                range(startPtNum + 1,len(coords))
        for endPtNum in endPtIndices:
            endPt = coords[endPtNum]
            dist = ((startPt[0] - endPt[0])**2 + 
                    (startPt[1] - endPt[1])**2)**(.5)
            if dist <= maxDist:
                thisNearest = dist
                break
            elif (dist < thisNearest) | (thisNearest == None):
                thisNearest = dist
                thisNearestID = endPtNum
        if (thisNearest > maxDist):
            farthest = [thisNearestID, startPtNum]
            maxDist = thisNearest
        startPtNum += 1

    print "Maximum distance to a nearest neighbor: " + str(maxDist)
    print "From " + str(farthest[0]) + " to " + str(farthest[1])
    return maxDist

if __name__ == "__main__":
    sys.exit(main())
