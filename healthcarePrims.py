#!/usr/bin/env python
"""
Uses Prim's Algorithm to build a network of a specified maximum number of nodes.
The network starts at the node closest to a given starting point. The starting
point is supplied as a single point in a separate shapefile.

Alex Zvoleff, aiz2101@columbia.edu
"""

import os
import sys
import getopt
import shutil
import time
from heapq import heappush, heappop

from modules import prims
from modules import fileRW
from modules import network

import arcgisscripting

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Error(Exception):
    def __init__(self, msg):
        self.msg = msg

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
        startPtShapefile = argv[2]
        outputShapefile = argv[3]
        nodeShapefile = argv[4]
        nodeShapefile1 = argv[5]
        lengthTable = argv[6]
        maxNodes = int(argv[7])

        gp.addMessage("Reading input nodes...")
        # Read nodes from input shapefile, creating network.Node instances.
        rows = gp.searchCursor(inputShapefile)
        row = rows.next()
        desc = gp.describe(inputShapefile)
        nodes = {}
        while row:
            feat = row.getValue(desc.ShapeFieldName)
            pnt = feat.getpart()
            # Read nodes.
            try:
                ptWeight = row.getValue("Weight")
            except:
                ptWeight = 1
            ptFID = row.getValue("FID")
            nodes[ptFID] = network.Node(ptFID, pnt.x, pnt.y, ptWeight)
            row = rows.next()
        del rows # ensure cursor closes

        outputPath = os.path.dirname(outputShapefile)
        if not os.path.exists(outputPath):
            os.mkdir(outputPath)
        distFile = outputPath + os.sep + "dists.txt"
        nodeFile = outputPath + os.sep + "nodes.txt"

        # Copy projection information to shpProj.prj file for later use
        rootDir, fc = os.path.split(inputShapefile)
        base, ext =  os.path.splitext(fc)
        projFile = rootDir + os.sep + base + ".prj"
        newProjFile = outputPath + os.sep + "shpProj.prj"
        try:
            if os.path.exists(newProjFile):
                os.remove(newProjFile)
            shutil.copy(projFile, newProjFile)
            projFile = newProjFile
        except WindowsError:
            gp.AddError("Could not find projection file")

        # Read starting point coordinates from the point in startPtShapefile
        rows = gp.searchCursor(startPtShapefile)
        desc = gp.describe(startPtShapefile)
        row = rows.next()
        feat = row.GetValue(desc.ShapeFieldName)
        FID = row.getValue("FID")
        point = feat.getPart(0)
        startNode = network.Node(FID, point.x, point.y, None)
        del rows
        
        firstNodeFID = nearestNode(startNode, nodes)

        gp.addMessage("Calculating tree...")
        tree = primsAlg(nodes, maxNodes, firstNodeFID)

        gp.addMessage("Writing output...")
        if os.path.exists(outputShapefile):
            gp.Delete_management(outputShapefile)
        fileRW.genShapefile(tree, projFile, outputShapefile)
        
        # Output a shapefile of the nodes in the network.
        if os.path.exists(nodeShapefile):
            gp.Delete_management(nodeShapefile)
        rootDir, fc = os.path.split(nodeShapefile)
        gp.CreateFeatureclass_management(rootDir, fc, "POINT")
        outDesc = gp.describe(nodeShapefile)
        shapefield = outDesc.ShapeFieldName
        origFIDField = "OrigFID"
        gp.addfield(nodeShapefile, origFIDField, "LONG") 
        lengthField = "Length"
        gp.addfield(nodeShapefile, lengthField, "FLOAT") 
        gp.deleteField_management(nodeShapefile,"Id")
        nodeRows = gp.insertCursor(nodeShapefile)
        point = gp.createobject("Point")
        for node in tree.getNodes():
            nodeFID = node.getID()
            sql = '"FID"=%i' %(nodeFID)
            rows = gp.searchCursor(inputShapefile, sql)
            oldRow = rows.next()
            row = nodeRows.newRow()
            row.setValue(origFIDField, nodeFID)
            row.setValue(lengthField, oldRow.getValue(lengthField))
            point.x = node.getX()
            point.y = node.getY()
            row.SetValue(shapefield, point)
            nodeRows.insertRow(row)
        del nodeRows
        del rows
        gp.defineprojection_management(nodeShapefile, projFile)
        if os.path.exists(nodeShapefile1):
            gp.Delete_management(nodeShapefile1)
        gp.CopyFeatures_management(nodeShapefile, nodeShapefile1)

        # Return a DBF with one field, the total length of the network.
        if os.path.exists(lengthTable):
            gp.Delete_management(lengthTable)
        rootDir, table = os.path.split(lengthTable)
        gp.CreateTable_management(rootDir, table)
        totalLenField = "totalLen"
        gp.addfield(lengthTable, totalLenField, "FLOAT") 
        gp.deleteField(lengthTable, "Field1")
        rows = gp.insertCursor(lengthTable)
        row = rows.newRow()
        row.setValue(totalLenField, tree.getTotalEdgeWeight())
        rows.insertRow(row)
        del rows
            
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "For help use --help"
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

def primsAlg(nodes, numNodes, firstNodeID = 0):
    'Prim\'s Algorithm for finding a minimum spanning tree'
    def addToHeap(heap, newSegs):
        'Adds new segments to the segHeap.'
        for seg in newSegs:
            heappush(heap, seg)
        return heap

    def calcSegs(startNode, nodes):
        """Calculates distance from startNode to each node in nodes, and returns a
        list of network.Seg instances"""
        segs = []
        for node in nodes:
            if node.getID() == startNode.getID():
                continue
            dist = ((startNode.getX() - node.getX())**2 + 
                    (startNode.getY() - node.getY())**2)**(.5)
            segs.append(network.Seg(None, startNode, node, dist))
        return segs

    tree = network.Network()
    segHeap = []

    # Find the shortest segment emanating from the node with the firstNodeID
    try:
        segs = calcSegs(nodes[firstNodeID], nodes.values())
    except KeyError:
        # If there are no segs from this node, there can be no tree.
        return tree

    leastWeight = None
    for seg in segs:
        if (seg.getWeight() < leastWeight) or (leastWeight == None):
            leastWeight = seg.getWeight()
            firstSeg = seg
    tree.addSeg(firstSeg)

    # Add the segs emanating from the first two endpoints to the heap
    for endNode in [firstSeg.getNode1(), firstSeg.getNode2()]:
        addToHeap(segHeap, calcSegs(endNode, nodes.values()))

    while tree.numNodes() < numNodes:
        try:
            seg = heappop(segHeap)
        except:
            # Tree is finished (not all nodes contained).
            break
        node1, node2 = seg.getNodes()
        node1InNet = tree.inNet(node1)
        node2InNet = tree.inNet(node2)
        # Add the seg if it's terminal node isn't already in the cluster.
        if (not node1InNet) or (not node2InNet):
            if not node1InNet:
                endNode = node1
            else:
                endNode = node2
            tree.addSeg(seg)
            # Add all emanating segs to the heap:
            addToHeap(segHeap, calcSegs(endNode, nodes.values()))
    return tree

if __name__ == "__main__":
    sys.exit(main())
