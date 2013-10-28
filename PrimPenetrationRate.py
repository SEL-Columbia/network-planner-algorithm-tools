#!/usr/bin/env python
"""

Expands an existing grid with different penetartion rates 

"""

import os
import sys
import getopt
import shutil
import pickle

import pylab
import numpy

import distStats
import distGen

import collections
import itertools

import calcTree

from modules import prims
from modules import fileRW
from modules import fileRWarcgis
from modules import network
from heapq import heappush, heappop
try:
    from osgeo import ogr
except ImportError:
    import ogr

try:
    from osgeo import osr
except ImportError:
    import osr

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Error(Exception):
    def __init__(self, msg):
        self.msg = msg


def readNodesFromShp(shapefile):
    'Reads nodes and node weights from a point shapefile.'
    ds = ogr.Open(shapefile)
    ptLayer = ds.GetLayer(0)
    nodes = {}
    feat = ptLayer.GetNextFeature()
    while feat is not None:
        nodeWeight = 1
        geomRef = feat.GetGeometryRef()
        x = geomRef.GetX()
        y = geomRef.GetY()
        FID = feat.GetFID()
        nodes[FID] = network.Node(FID, x, y, nodeWeight)
        feat = ptLayer.GetNextFeature()
    ds.Destroy()

    return nodes

def generateDistFromNodeDict(nodeDict,inputFile, outputPath, searchRadius):
    searchRadius = fileRW.convertToShpUnits(inputFile, searchRadius)
    rootDir, fc = os.path.split(inputFile)
    file, ext = os.path.splitext(fc)

    if not os.path.exists(outputPath):
        try:
            os.mkdir(outputPath)
        except:
            print "ERROR: could not create new directory", outputPath

    print "Calculating distances..."
    dists = []
    nodeCopy = nodeDict.copy()
    for startNode in nodeDict.values():
        del nodeCopy[startNode.getID()]
        for endNode in nodeCopy.values():
            dist = ((startNode.getX() - endNode.getX())**2 + 
                    (startNode.getY() - endNode.getY())**2)**(.5)
            if dist < searchRadius:
                dists.append([startNode.getID(), endNode.getID(), dist])

    print "Outputting nodes..."
    nodeFile = outputPath + os.sep + "nodes.txt"
    fileRW.writeNodesToTxt(nodeDict, nodeFile)

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

def primsAlgforLV(nodes,tree, numNodes, nodeDict, numberOfNodesCovered, excludedNodes = []): 
    'Prim\'s Algorithm for finding a minimum spanning tree'
    ### nodeDict=segmentsByNodeID

    segHeap = []
    for endNode in tree.getNodes():
        addToHeap(segHeap, nodeDict[endNode.getID()])
    
    while tree.numNodes() < (numberOfNodesCovered):
        try:
            # Get best segment from heap
            seg = heappop(segHeap)
        except:
            # Tree is finished (not all nodes contained).
            break
        node1, node2 = seg.getNodes() 
        node1InNet = tree.inNet(node1)
        node2InNet = tree.inNet(node2)
        # Add the seg if it's terminal node isn't already in the cluster.
        if (not node1InNet) or (not node2InNet):
        #if ((not node1InNet) and node2InNet) or ((not node2InNet) and node1InNet):
            
            if not node1InNet:
                endNode = node1
            else:
                endNode = node2
            tree.addSeg(seg)
            # Add all emanating segs to the heap:
            # nodeDict returns all segments coming out from the endNode
            # endNode is the node that is outside of the tree
            addToHeap(segHeap, nodeDict[endNode.getID()])
            # And we are sure that everything in the heap is adjacent to the tree because
            # we only add the adjacent segments in the first place using nodeDict
    return tree

def addToHeap(heap, newSegs):
    'Adds new segments to the segHeap.'
    for seg in newSegs:
        heappush(heap, seg)
    return heap

def buildAssocDict(segments, excludedNodes = []):
    'Builds dictionary with nodeID key where values are segs from/to that node'
    segList = {}
    for seg in segments:
        node1, node2 = seg.getNodes()
        if (node1 in excludedNodes) or (node2 in excludedNodes):
            continue
        else:
            for nodeID in [node1.getID(), node2.getID()]:
                if segList.has_key(nodeID):
                    segList[nodeID].append(seg)
                else:
                    segList[nodeID] = [seg]
    return segList


                                                    

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.GetoptError, msg:
            raise Usage(msg) 

        try:
            #nodeFile = argv[1] # Nodes Shape File
            #nodeFile=r"C:\Users\Selin\Documents\COLUMBIA\Project\for World Bank Prensentation\RuhiiraForVillageNetwork\AllPointsIncludingFakes.shp"
            nodeFile=r"C:\Users\Selin\Documents\COLUMBIA\Project\for World Bank Prensentation\RuhiiraForVillageNetwork\MainGridExpansion\AllNodesMerged_Updated.shp"
            #OutputDir = argv[2] # Directory for texts and final network
            OutputDir=r"C:\Users\Selin\Documents\COLUMBIA\Project\for World Bank Prensentation\RuhiiraForVillageNetwork\MainGridExpansion\75"
            #existingTree=argv[3]
            #existingTree=r"C:\Users\Selin\Documents\COLUMBIA\Project\for World Bank Prensentation\RuhiiraForVillageNetwork\Backbone1\network.shp"
            existingTree=r"C:\Users\Selin\Documents\COLUMBIA\Project\for World Bank Prensentation\RuhiiraForVillageNetwork\MainGridExpansion\backBone\network.shp"
            searchRadius = "1500 meters"
            #numNodes=8912

        except IndexError:
            raise Error("Not enough arguments provided to script.")
        
        nodes=readNodesFromShp(nodeFile)
        numNodes=len(nodes)
        net=fileRWarcgis.readNetFromShpUsingNodes(existingTree, nodes)
        backBoneNodes=net.numNodes()
        print backBoneNodes, numNodes

                
        generateDistFromNodeDict(nodes,nodeFile,OutputDir, searchRadius) ### bu functiona dikkat!!!proj file!!!
        
        shapePath = OutputDir + os.sep + "ExpandedNetwork.shp"
        distFileLV = OutputDir + os.sep + "dists.txt"
        nodeFileLV = OutputDir + os.sep + "nodes.txt"
        projFile = OutputDir + os.sep + "shpProj.prj"
        testSegmentsForLV=fileRW.readSegsFromTxt(distFileLV, nodes)
        segmentsByNodeID=buildAssocDict(testSegmentsForLV, excludedNodes=[])
        resultsFile = OutputDir + os.sep + "results.txt"
        ofile = open(resultsFile, "w")
        ofile.write("PenetrationRate NodesConnected Total Lenght\n")
        for i in range(10,11):
            penetrationRate=0.05*i
            finalTree=primsAlgforLV(nodes,net,numNodes, segmentsByNodeID,(backBoneNodes+penetrationRate*(numNodes-backBoneNodes)),excludedNodes=[])
            totalCost=finalTree.getTotalEdgeWeight()
            nodesConnected=finalTree.numNodes()
            ofile.write("%(penetrationRate)f %(nodesConnected)i %(totalCost)f \n" %vars())  
        
        ofile.write("END\n")
        ofile.close()
        
        print finalTree.getTotalEdgeWeight()   
        fileRW.genShapefile(finalTree, projFile, shapePath)
        
    except Usage, err:
            
        print >>sys.stderr, err.msg
        return 2

    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

if __name__ == "__main__":
    sys.exit(main())


