#!/usr/bin/env python
"""
Uses prim's algorithm to build a series of subnets of maximum length (in nodes)
specified by an input parameter. The subnetIDs and order (within each subnet)
are written to the input shapefile. The networks are calculated iteratively,
with the first point chosen by choosing the point closest to a given starting
point. After each subnet is formed, those nodes are removed from the set of
input points, and the next subnet calculated.

Alex Zvoleff, aiz2101@columbia.edu
"""
#TODO: Not finished. Complete the coding of finding the starting node.
#TODO: Rewrite fileRW.genShapeFile to handle a list of trees

import os
import sys
import getopt

import distGen
from heapq import heappush, heappop

from modules import fileRW
from modules import network

import arcgisscripting

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Error(Exception):
    def __init__(self, msg):
        self.msg = msg

class distance():
    def __init__(self, dist, node):
        self._dist = dist
        self._node = node

    def __lt__(self, other):
        if self._dist < other._dist:
            return True
        else:
            return False

    def __eq__(self, other):
        if self._node == other._node:
            return True
        else:
            return False

    def __gt__(self, other):
        if self._dist > other._dist:
            return True
        else:
            return False

    def getNode(self):
        return self._node

try:
    gp = arcgisscripting.create()
except:
    sys.exit("ERROR: Could not load geoprocessor. Check license server.")

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

        inputShapefile = argv[1]
        distancesDir = argv[2]
        startPtShapefile = argv[3]
        maxNodes = float(argv[4])
        subnetIDFieldname = argv[5]
        orderFieldname = argv[6]

        distFile = distancesDir + os.sep + "dists.txt"
        nodeFile = distancesDir + os.sep + "nodes.txt"
        nodes = fileRW.readNodesFromTxt(nodeFile)
        segs = fileRW.readSegsFromTxt(distFile, nodes)
        
        gp.addMessage("Reading starting point...")
        rows = gp.searchCursor(startPtShapefile)
        desc = gp.describe(startPtShapefile)
        row = rows.next()
        feat = row.GetValue(desc.ShapeFieldName)
        FID = row.getValue("FID")
        point = feat.getPart(0)
        startNode = network.Node(FID, point.x, point.y, None)
        del rows

        nodeDict = buildAssocDict(segs)

        gp.addMessage("Calculating subnets...")
        nodeOrder = {}
        subNet = {}
        distList = makeDistList(startNode, nodes)
        n = 0
        while len(segs) > 0:
            print n
            # Choose first node
            firstNodeID = distList[0].getNode()
            tree = primsAlg(segs, maxNodes, firstNodeID, nodeDict)
            m = 0
            for node in tree.getNodes():
                nodeID = node.getID()
                print nodeID
                nodeOrder[nodeID] = m
                subNet[nodeID] = n
                distList.remove(distance(None,nodeID))
                for seg in nodeDict[node.getID()]:
                    for node in [seg.getNode1(), seg.getNode2()]:
                        nodeDict[node.getID()].remove(seg)
                m += 1
            n += 1
            if n == 5:
                break

        gp.addMessage("Writing data to shapefile...")
        # Copy projection information to shpProj.prj file for later use
        rootDir, fc = os.path.split(inputShapefile)
        base, ext =  os.path.splitext(fc)
        projFile = rootDir + os.sep + base + ".prj"
        outputShapefile = rootDir + os.sep + base + "nets" + ".shp"
        fileRW.genShapefile(tree, projFile, outputShapefile)
        #fileRW.writeFieldToShp(inputShapefile, nodeOrder, orderFieldname)
        #fileRW.writeFieldToShp(inputShapefile, subNet, subnetIDFieldname)

    
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "For help use --help"
        return 2

    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

def makeDistList(startNode, nodes):
    """Makes a list of distances to simplify finding the starting point at each
    iteration"""
    distList = []
    for node in nodes.values():
        dist = ((startNode.getX() - node.getX())**2 + 
                (startNode.getY() - node.getY())**2)**(.5)
        distList.append(distance(dist, node.getID()))
    distList.sort()
    return distList

def primsAlg(segments, numNodes, firstNodeID, nodeDict):
    'Prim\'s Algorithm for finding a minimum spanning tree'
    tree = network.Network()
    segHeap = []

    # Find the shortest segment emanating from the node with the firstNodeID
    try:
        segs = nodeDict[firstNodeID]
    except KeyError:
        #When the tree is only this node... no emanating segs
        return tree

    leastWeight = None
    for seg in segs:
        if (seg.getWeight() < leastWeight) or (leastWeight == None):
            leastWeight = seg.getWeight()
            firstSeg = seg
    tree.addSeg(firstSeg)

    # Add the segs emanating from the first two endpoints to the heap
    for endNode in [firstSeg.getNode1(), firstSeg.getNode2()]:
        newSegs = nodeDict[endNode.getID()]
        addToHeap(segHeap,newSegs)

    while tree.numNodes() < numNodes:
        try:
            seg = heappop(segHeap)
        except:
            # Tree is finished (not all nodes contained).
            break

        node1 = seg.getNode1()
        node2 = seg.getNode2()
        node1InNet = tree.inNet(node1)
        node2InNet = tree.inNet(node2)

        # Add the seg if it's terminal node isn't already in the cluster.
        if (node1InNet != node2InNet):
            if not node1InNet:
                endNode = node1
            else:
                endNode = node2
            tree.addSeg(seg)
            # Add all emanating segs to the heap:
            newSegs = nodeDict[endNode.getID()]
            addToHeap(segHeap,newSegs)
    return tree

def addToHeap(heap,newSegs):
    'Adds new segments to the segHeap.'
    for seg in newSegs:
        heappush(heap,seg)
    return heap

def buildAssocDict(segments):
    'Builds dictionary with nodeID key where values are segs from/to that node'
    segList = {}
    for seg in segments:
        node1 = seg.getNode1()
        node2 = seg.getNode2()
        for nodeID in [node1.getID(), node2.getID()]:
            if segList.has_key(nodeID):
                segList[nodeID].append(seg)
            else:
                segList[nodeID] = [seg]
    return segList

if __name__ == "__main__":
    sys.exit(main())
