#!/usr/bin/env python
"""
Runs Kruskal, breask before last (ClusterNumber-1) segment.

"""

import os
import sys
import getopt
import shutil
import time
import pickle

import distStats
import distGen

import collections
import itertools

import calcTree

from modules import prims
from modules import fileRW
from modules import network
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
    nodesArray=[]
    feat = ptLayer.GetNextFeature()
    while feat is not None:
        #try:
         #   weightField = feat.GetFieldIndex("Weight")
          #  nodeWeight = feat.GetField(weightField)
        #except:
        nodeWeight = 1
        geomRef = feat.GetGeometryRef()
        x = geomRef.GetX()
        y = geomRef.GetY()
        FID = feat.GetFID()
        nodes[FID] = network.Node(FID, x, y, nodeWeight)
        nodesArray.append((x,y))
        feat = ptLayer.GetNextFeature()
    ds.Destroy()
    return nodes, nodesArray


def meanCluster(nodeArray,centerNumber):
    results = Pycluster.kcluster(nodeArray,nclusters=centerNumber,npass=50,method='m')
    assignments=results[0]
    xs = [node[0] for node in nodeArray]
    ys = [node[1] for node in nodeArray]
    clusterByClusterID=collections.defaultdict(list)
    #clusterCenterByClusterID=collections.defaultdict(list)
    clusterCenterArray=[]
    for x, y, clusterID in itertools.izip(xs,ys,assignments):
        clusterByClusterID[clusterID].append((x,y))

    for clusterID in clusterByClusterID:
        mean=numpy.mean(clusterByClusterID[clusterID],axis=0)
        clusterCenterArray.append(mean)
    return clusterCenterArray, clusterByClusterID

def kruskalForCluster(segments, numNodes,centerNumber, maxSegWeight=None):
    "Kruskal\'s algorithm for finding clusters"
    segments.sort(key=lambda obj:obj.getWeight())
    tree = network.Network()
    
    for segment in segments:
        node1 = segment.getNode1()
        node2 = segment.getNode2()
        node1InNet = tree.inNet(node1)
        node2InNet = tree.inNet(node2)
        
        if (not node1InNet and not node2InNet) or (node1InNet != node2InNet): 
            tree.addSeg(segment)
        else:
             if node1InNet and node2InNet and \
                    (tree.getNetID(node1) != tree.getNetID(node2)):
                        tree.addSeg(segment)
        if tree.numNodes() > numNodes-centerNumber+1:
            break
        if (maxSegWeight is not None) and (segment.getWeight() > maxSegWeight):
            break
    return tree




def runAlgorithm(algorithm, inputDir, outputDir,centerNumber, cutoff=None):
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
    
    print "\tCalculating tree..."
    
    tree = kruskalForCluster(testSegments, numNodes,centerNumber)

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


                                                    

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.GetoptError, msg:
            raise Usage(msg) 

        try:
            inputDir = argv[1] # Directory for dist and nodes text files.
            inputDir="C:\Documents and Settings\Selin\Desktop\Mbola\MbolaPts.shp"
            outputDir=argv[2] #Directory for output shapefile
            outputDir="C:\Documents and Settings\Selin\Desktop\Nov21Cluster"
            centerNumber=10
            algorithm="Kruskal's Algorithm" # Kruskal for Cluster"
            numNodes=1174  #### Make this automated
        except IndexError:
            raise Error("Not enough arguments provided to script.")

    runAlgorithm(algorithm, inputDir, outputDir,centerNumber, cutoff)
    except Usage, err:
            
        print >>sys.stderr, err.msg
        return 2

    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

if __name__ == "__main__":
    sys.exit(main())


