#!/usr/bin/env python
"""
Trims a precalculated tree by removing the longest segment until the
specified number of nodes is reached.

Alex Zvoleff, aiz2101@columbia.edu
"""
#TODO: Adapt to work with new network class, and test.

import os
import sys
import getopt
import heapq

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

        shapefile = argv[1]
        fieldName = argv[2]

        orderCut = treeTrim(shapefile)
        fileRW.writeFieldToShp(shapefile, orderCut, fieldName)

    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "For help use --help"
        return 2

    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

def treeTrim(shapefile):
    """Progressively removes the longest segment from a tree, storing the order
    in which the segments are removed"""
    tree = fileRW.readNetFromShp(shapefile)
    segs = tree.getEdges()
    nodes = tree.getNodes()

    numNodes = len(nodes)
    nodeDict = buildAssocDict(segs)
    terminalNodes = []
    for node in nodeDict.iterkeys():
        if len(nodeDict[node]) == 1:
            terminalNodes.append(nodeDict[node][0])
    terminalHeap = []
    addToHeap(terminalHeap,terminalNodes)
    orderCut = {} #track order segments are removed
    n = 0
    while len(nodeDict) > 0:
        cutSeg = heapq.nlargest(1,terminalHeap)[0]
        terminalHeap.remove(cutSeg)
        orderCut[cutSeg.getID()] = n
        endPts = [cutSeg.getNode1(), cutSeg.getNode2()]
        for endPt in endPts:
            if len(nodeDict[endPt]) == 1:
                del nodeDict[endPt]
            else:
                endPtIndex = nodeDict[endPt].index(cutSeg)
                del nodeDict[endPt][endPtIndex]
                if len(nodeDict[endPt]) == 1:
                    addToHeap(terminalHeap,nodeDict[endPt])
        n += 1
    return orderCut

def addToHeap(heap, newSegs):
    'Adds new segments to the terminalHeap.'
    for seg in newSegs:
            heapq.heappush(heap, seg)
    return heap

def buildAssocDict(segments):
    'Builds dictionary with nodeID key where values are segs from/to that node'
    segList = {}
    for seg in segments:
        node1 = seg.getNode1()
        node2 = seg.getNode2()
        for node in [node1, node2]:
            if segList.has_key(node):
                segList[node].append(seg)
            else:
                segList[node] = [seg]
    return segList

if __name__ == "__main__":
    sys.exit(main())
