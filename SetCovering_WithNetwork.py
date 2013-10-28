#!/usr/bin/env python
"""
A Set Covering Heuristic Algorithm
Decides number and the locations of the facilities.
Ayse Selin Kocaman
ask2170@columbia.edu
"""

import os
import sys
import getopt
import time
import copy
import CMST_dfs
import gc
import collections
import scipy
import csv
#import pylab
import numpy
#import batchPrimsforTransformers
from heapq import heappush, heappop
from osgeo import ogr
from modules import network
from modules import fileRW
from modules import prims
from collections import defaultdict

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Error(Exception):
    def __init__(self, msg):
        self.msg = msg   
     
def generateNodeDictFromShp(shapeFile,outputPath): 
    'Reads nodes and node weights from a point shapefile.'
    rootDir, fc = os.path.split(shapeFile)
    file, ext = os.path.splitext(fc)

    if not os.path.exists(outputPath):
        try:
            os.mkdir(outputPath)
        except:
            print "ERROR: could not create new directory", outputPath
    ds = ogr.Open(shapeFile)
    ptLayer = ds.GetLayer(0)
    
    #nodesByClusterID=collections.defaultdict(list)
    #clustersByNodeID={}
    nodes={}
       
    feat = ptLayer.GetNextFeature()
    while feat is not None:
        nodeWeight = 1
        geomRef = feat.GetGeometryRef()
        x = geomRef.GetX()
        y = geomRef.GetY()
        FID = feat.GetFID()
        nodes[FID] = network.Node(FID, x, y, nodeWeight) #Households
        feat = ptLayer.GetNextFeature()
    ds.Destroy()
    return nodes

def generateClusterDicts(nodes,coverDist): 
    nodesByClusterID=collections.defaultdict(list)
    clustersByNodeID=collections.defaultdict(list)
   
    for startNode in nodes.values():
        for endNode in nodes.values():
            dist = ((startNode.getX() - endNode.getX())**2 + 
                    (startNode.getY() - endNode.getY())**2)**(.5)
            if dist < coverDist:
                nodesByClusterID[startNode.getID()].append(endNode.getID())
                clustersByNodeID[endNode.getID()].append(startNode.getID())
                
    return nodesByClusterID,clustersByNodeID


def maxInClusterDist(centerNode,nodesByClusterID): #Returns maxDist within the cluster
    maxdist=0
    for node in nodesByClusterID[centerNode.getID()]: #uses the fact that centerID and ClusterID are same
        dist=((centerNode.getX()-node.getX())**2+
                (centerNode.getY()-node.getY())**2)**(.5)
        if dist>=maxdist:
            maxdist=dist
    return maxdist



    
def totalInClusterCost(nodesByClusterID,centers):
    totalCost=0
    for centerID in centers.keys():
        for node in nodesByClusterID[centerID]:
            totalCost+=((node.getX()-centers[centerID].getX())**2+
                        (node.getY()-centers[centerID].getY())**2)**(.5)
    return totalCost
    
def addLVSeg(tree,centers,nodesByClusterID):#single points line from the root
    
    SegID=1000000
    meanSeg={}
    for centerID in centers.keys():
        tree._nodesByNetID[centerID]=[]
        tree._network[centerID]=[]
        sumLen=0
        count=0
        for node in nodesByClusterID[centerID]:
            length=((node.getX()-centers[centerID].getX())**2+
                        (node.getY()-centers[centerID].getY())**2)**(.5)
            newSeg=network.Seg(SegID, node, centers[centerID], length)
            tree._netIDByNode[node] = centerID
            tree._nodesByNetID[centerID].append(node)
            tree._network[centerID].append(newSeg)
            sumLen+=length
            count+=1
        mean= sumLen/count  
        meanSeg[centerID]=mean
    #print tree.getEdges()
    return tree,meanSeg

'''def writeCenterSizeToText(statsFile, Dict):
    # import csv
    # csvWriter = csv.writer(open(statsFile, 'wb'))
    # for key in Dict:
    #    csvWriter.writerow([Dict[key].getX(), Dict[key].getY(), Dict[key].getWeight()])

   
    outFile = open(statsFile,"w")
    for key in Dict.keys():
            size=Dict[key].getWeight()
	    x=Dict[key].getX()
	    y=Dict[key].getY()
            outFile.write("%(size)i \n" %vars())
	    outFile.write("%(x)f \n" %vars())
	    outFile.write("%(y)f \n" %vars())
    outFile.close()
    return 0'''


def findTheBiggestCluster(nodesByClusterID):
    maxID=0
    maxLen=0
    for nodeID in nodesByClusterID.keys():
        l=len(nodesByClusterID[nodeID])
        if l>maxLen:
            maxLen=l
            maxID=nodeID
    return maxID

def updateDicts(nodesByClusterID,maxID,clusters,nodes):
   
    copy=[]
    
    for nodeID in nodesByClusterID[maxID]: #burayi bir kontrol et bakalim
        copy.append(nodeID)
        clusters[maxID].append(nodes[nodeID])
        del nodesByClusterID[nodeID]
   
    for nodeID in nodesByClusterID.keys():
        copy2=[]
        for n in nodesByClusterID[nodeID]:
            copy2.append(n)
        for node in copy2:
            for ID in copy:
                if node==ID:
                    nodesByClusterID[nodeID].remove(node)
        del copy2           
    del copy
    return nodesByClusterID,clusters

def kruskalsAlg(segments,nodes):
    'Kruskal\'s algorithm for finding a minimum spanning tree'
    segments.sort(key=lambda obj:obj.getWeight())
    tree = network.Network()
    numNodes=len(nodes)
    
    for segment in segments:
        node1 = segment.getNode1()
        node2 = segment.getNode2()
        node1InNet = tree.inNet(node1)
        node2InNet = tree.inNet(node2)
        
        if (not node1InNet and not node2InNet) or (node1InNet != node2InNet): 
            tree.addSeg(segment)
           
        else:
             if node1InNet and node2InNet and (tree.getNetID(node1) != tree.getNetID(node2)):
                tree.addSeg(segment)
                
        if tree.numNodes() > numNodes:
            break

    return tree,segments

def generateSegments(centers,searchRadius): #centers=centerDict
    segments=[]
    #nodeCopy = centers.copy()
    nodeCopy=copy.deepcopy(centers)
    segID=0
    for startNode in centers.values():
        del nodeCopy[startNode.getID()]
        for endNode in nodeCopy.values():
            dist = ((startNode.getX() - endNode.getX())**2 + 
                    (startNode.getY() - endNode.getY())**2)**(.5)
            if dist < searchRadius:
                segments.append(network.Seg(segID, startNode, endNode, dist))
                segID+=1
    return segments

def generateRootNode(nodesByClusterID, ID, nodes):
    sumX=0
    sumY=0
    numNodes=0
    for nodeID in nodesByClusterID[ID]:
        sumX=sumX+nodes[nodeID].getX()
        sumY=sumY+nodes[nodeID].getY()
        numNodes+=1
        

    centerX=sumX/numNodes
    centerY=sumY/numNodes
    newNode=network.Node(ID, centerX, centerY, 1) # key can be changed! karismasin simdilik        
    return newNode

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.GetoptError, msg:
            raise Usage(msg)
        try:
            #inputShapeFile = argv[1]
            
            inputShapeFile=r"/home/selin/Desktop/windows-desktop/real/Ikaram/IkaramPts.shp"
            
            #outputDir = argv[2]
            outputDir=r"/home/selin/Desktop/windows-desktop/real/Ikaram/IkaramPts"
            SR=100000000
            coverDist = 500
            maxLVLenghtInCluster=600
        except IndexError:
            raise Error("Not enough arguments provided to script.")

        startTime = time.time()

        clusters=collections.defaultdict(list)
        print "Generating Dictionaries"
        nodes=generateNodeDictFromShp(inputShapeFile,outputDir)
       
        nodesByClusterID,clustersByNodeID=generateClusterDicts(nodes,coverDist)
        
        statsFile= outputDir + os.sep + "Centers_SetCover.csv"
        csvWriter = csv.writer(open(statsFile, 'wb'))
        outFile = open(statsFile,"w")
        centers={}
        while nodesByClusterID:
            maxID=findTheBiggestCluster(nodesByClusterID)
            centers[maxID]=nodes[maxID]
            #centers[maxID]=generateRootNode(nodesByClusterID, maxID,nodes) #eger center of mass istersem
            size=len(nodesByClusterID[maxID])
            nodesByClusterID,clusters=updateDicts(nodesByClusterID,maxID,clusters,nodes)
            
            x=nodes[maxID].getX()
            y=nodes[maxID].getY()
            csvWriter.writerow([x, y, size])
        outFile.close()
        print "Num Clusters", len(clusters)
        print "Total Running Time:",time.time()-startTime

        '''
        try:
            netID=tree.getNetID(centers.values()[0])
        except:
            netID=0
            tree._nodesByNetID[0]=[]
            tree._network[netID]=[]
        '''


        
        
        #roots={}
        #roots=generateRootNode(nodesByClusterID)
        
        tree=network.Network()

        print "Generating LV network"
        
        netID=0
        tree._nodesByNetID[0]=[]
        tree._network[netID]=[]

       
            
        
        for ID in clusters.keys():
            nodesByNodeID={}
            segments,lvCost=CMST_dfs.CMST(clusters[ID],maxLVLenghtInCluster,centers[ID])
           
            for segment in segments.values():
                node1=segment.getNode1()
                node2=segment.getNode2()
                if not nodesByNodeID.has_key(node1.getID()):
                    nodesByNodeID[node1.getID()]=node1

                if not nodesByNodeID.has_key(node2.getID()):
                    nodesByNodeID[node2.getID()]=node2

            for node in nodesByNodeID.values():                 
                tree._netIDByNode[node] = netID
                tree._nodesByNetID[netID].append(node)

            for segment in segments.values():
                tree._network[netID].append(segment)

        print "Generating MV network"
        

        LVLen=tree.getTotalEdgeWeight()
        print "LVLen", LVLen
        segments=generateSegments(centers,SR)
        
        MVtree,segments=kruskalsAlg(segments,centers)
        MVLen=MVtree.getTotalEdgeWeight()
        print "MVLen", MVLen
        fileRW.genShapefile(MVtree, outputDir + ".prj", outputDir + os.sep + "MV.shp")
        fileRW.genShapefile(tree, outputDir + ".prj", outputDir + os.sep + "LV.shp")
       
        
    except Usage, err:
        print >>sys.stderr, err.msg

        return 2
    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

if __name__ == "__main__":
    sys.exit(main())

