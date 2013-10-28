#!/usr/bin/env python
"""
A heuristic for facility location algorithm
Ayse Selin Kocaman
ask2170@columbia.edu
"""

import os
import sys
import getopt
import time
import copy
#import CMST_dfs
import gc
#import pylab
import numpy
#import batchPrimsforTransformers
import collections
from osgeo import ogr
from modules import network
from modules import fileRW
from modules import prims

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Error(Exception):
    def __init__(self, msg):
        self.msg = msg
   
def generateDictsFromShp(shapeFile,outputPath): 
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
    #clusters = {}
    nodesByClusterID=collections.defaultdict(list)
    clusterByNode={}
    nodes={}
    centers={}
    LVCostDict={}
    #nodesArray=[]
    feat = ptLayer.GetNextFeature()
    while feat is not None:
        nodeWeight = 1
        geomRef = feat.GetGeometryRef()
        x = geomRef.GetX()
        y = geomRef.GetY()
        FID = feat.GetFID()
        nodes[FID] = network.Node(FID, x, y, nodeWeight) #tum nodelar
        centers[FID]=network.Node(FID, x, y, nodeWeight) #tum centerlar node olarak tutulur
        #clusters[FID]=Cluster([nodes[FID]],FID,x,y,nodeWeight) #tum clusterlar
        
        #clusterByNode[nodes[FID]]=Cluster([nodes[FID]],FID,x,y,nodeWeight) # node'u ver clusteri dondur
        clusterByNode[nodes[FID]]=FID #clusterin kendini degil sadece IDsini tutucak
        #nodesByClusterID[FID]=clusterByNode[nodes[FID]].getNodes() #clusterID ver clusterin nodelarini dondur
        nodesByClusterID[FID].append(nodes[FID])
        LVCostDict[FID]=0
        #nodesArray.append((x,y))
        feat = ptLayer.GetNextFeature()
    ds.Destroy()
    return nodesByClusterID,clusterByNode,nodes,centers,LVCostDict

def generateSegments(centers,searchRadius): #centers=centerDict
    #segments = {}
    segments=[]
    nodeCopy = centers.copy()
    segID=0
    for startNode in centers.values():
        del nodeCopy[startNode.getID()]
        for endNode in nodeCopy.values():
            dist = ((startNode.getX() - endNode.getX())**2 + 
                    (startNode.getY() - endNode.getY())**2)**(.5)
            if dist < searchRadius:
                
                segments.append(network.Seg(segID, startNode, endNode, dist))
                #segments[(startNode.getID(), endNode.getID())]=network.Seg(segID,startNode,endNode,dist)
                segID+=1
    return segments



def maxInClusterDist(centerNode,nodesByClusterID): #Returns dist array within the cluster
    maxdist=0
    for node in nodesByClusterID[centerNode.getID()]: #uses the fact that centerID and ClusterID are same
        dist=((centerNode.getX()-node.getX())**2+
                (centerNode.getY()-node.getY())**2)**(.5)
        if dist>=maxdist:
            maxdist=dist
    return maxdist

def maxTempInClusterDist(segment,ClusterByNode,nodesByClusterID):
    maxDist=0

    tempCenter1,tempCenter2=segment.getNodes()
       
    tempCenterX=(tempCenter1.getWeight()*tempCenter1.getX()
        +tempCenter2.getWeight()*tempCenter2.getX())/(tempCenter2.getWeight()+tempCenter1.getWeight())
    tempCenterY=(tempCenter1.getWeight()*tempCenter1.getY()
        +tempCenter2.getWeight()*tempCenter2.getY())/(tempCenter2.getWeight()+tempCenter1.getWeight())
    
    for node in nodesByClusterID[ClusterByNode[segment.getNode1()]]:
        dist=((tempCenterX-node.getX())**2+(tempCenterY-node.getY())**2)**(.5)
        if dist>=maxDist:
            maxDist=dist

    
    for node in nodesByClusterID[ClusterByNode[segment.getNode2()]]:
        dist=((tempCenterX-node.getX())**2+(tempCenterY-node.getY())**2)**(.5)
        if dist>=maxDist:
            maxDist=dist

    return maxDist,tempCenterX,tempCenterY

def totalInClusterCost(nodesByClusterID,centers):
    totalCost=0
    for centerID in centers.keys():
        for node in nodesByClusterID[centerID]:
            totalCost+=((node.getX()-centers[centerID].getX())**2+
                        (node.getY()-centers[centerID].getY())**2)**(.5)
    return totalCost
    
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
             if node1InNet and node2InNet and \
                    (tree.getNetID(node1) != tree.getNetID(node2)):
                        tree.addSeg(segment)
        if tree.numNodes() > numNodes:
            break
        
    return tree,segments
def totalInClusterDistance(nodesByClusterID,centers):
    totalDistance=0
    for centerID in centers.keys():
        for node in nodesByClusterID[centerID]:
            totalDistance+=((node.getX()-centers[centerID].getX())**2+
                        (node.getY()-centers[centerID].getY())**2)**(.5)
    return totalDistance


   
def run(centers,nodesByClusterID,clusterByNode,sr): 
    segments=generateSegments(centers,sr)
    minTree,segments=kruskalsAlg(segments,centers) #tree1 is minTree at the beginning.
    minTotalCost=len(centers)*10000
    minNodesByClusterID=copy.deepcopy(nodesByClusterID)
    minCenters=copy.deepcopy(centers)
    
    #fix tempCenterX , tempCenterY and maxDist of the first merge before going into while
    if segments[0].getWeight()<=200:
        maxDist=0 # can be anything less than 500
    else:
        maxDist=101 # can be anything greater than 500
        print "NO CLUSTER POSSIBLE"

    tempCenter1,tempCenter2=segments[0].getNodes()
       
    tempCenterX=(tempCenter1.getWeight()*tempCenter1.getX()
        +tempCenter2.getWeight()*tempCenter2.getX())/(tempCenter2.getWeight()+tempCenter1.getWeight())
    tempCenterY=(tempCenter1.getWeight()*tempCenter1.getY()
        +tempCenter2.getWeight()*tempCenter2.getY())/(tempCenter2.getWeight()+tempCenter1.getWeight())
    
    while(maxDist<=100):    

        segment=segments[0]
        
        # merge 
        center1,center2=segment.getNodes()

        weight=center2.getWeight()+center1.getWeight()
        baseClusterID=min(clusterByNode[center1],clusterByNode[center2])
        mergingClusterID=max(clusterByNode[center1],clusterByNode[center2])
   
        nodesByClusterID[baseClusterID].extend(nodesByClusterID.pop(mergingClusterID))
    
        centers[baseClusterID].setXY(tempCenterX,tempCenterY)
        centers[baseClusterID].setWeight(weight)
    
        del centers[mergingClusterID] #yeni center ve cluster eklemek lazim
   
        for node in nodesByClusterID[baseClusterID]:
            clusterByNode[node]=baseClusterID

        segments=generateSegments(centers,sr) # generate segments for new graph
        newTree,segments=kruskalsAlg(segments,centers)#returns sorted segments
       
        TotalTransformerCost=len(centers)*10000
       
        LVCost=totalInClusterDistance(nodesByClusterID,centers)*2
        newTotalCost=TotalTransformerCost+LVCost

        if(newTotalCost<=minTotalCost):
            minNodesByClusterID=copy.deepcopy(nodesByClusterID)
            minTree=copy.deepcopy(newTree)
            minCenters=copy.deepcopy(centers)
            minTotalCost=newTotalCost

        # Calculate maxDist below for next graph and continue if it is less than 500

        try:
            segment=segments[0]
        except:
            break   
        maxDist,tempCenterX,tempCenterY=maxTempInClusterDist(segment,clusterByNode,nodesByClusterID)
           
    return minTotalCost,minTree,minCenters,minNodesByClusterID


def addLVSeg(tree,centers,nodesByClusterID):#single points line from the root
    
    SegID=1000000
    
    for centerID in centers.keys():
        tree._nodesByNetID[centerID]=[]
        tree._network[centerID]=[]
       
        for node in nodesByClusterID[centerID]:
            length=((node.getX()-centers[centerID].getX())**2+
                        (node.getY()-centers[centerID].getY())**2)**(.5)
            newSeg=network.Seg(SegID, node, centers[centerID], length)
            tree._netIDByNode[node] = centerID
            tree._nodesByNetID[centerID].append(node)
            tree._network[centerID].append(newSeg)
    #print tree.getEdges()
    return tree


def writeCentersToText(statsFile1,statsFile2, Dict):
    'Writes Center dict to a text file for future analysis'
    outFile1 = open(statsFile1,"w")
    outFile2 = open(statsFile2,"w")
    for key in Dict.keys():
        id = Dict[key].getID()
        x, y = Dict[key].getX(), Dict[key].getY()
        weight = Dict[key].getWeight()
        outFile1.write("%(weight)i\n" %vars())
        outFile2.write("%(id)i %(x).11e %(y).11e %(weight)f\n" %vars())
    outFile1.close()
    outFile2.close()
    return 0
   
def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.GetoptError, msg:
            raise Usage(msg)

        #if len(args) == 0:
            #raise Usage("No arguments specified.")
        try:
            #inputShapeFile = argv[1]
            #inputShapeFile=r"C:\Documents and Settings\Selin\Desktop\Mbola\MbolaPts.shp"
            inputShapeFile=r"C:\Documents and Settings\Selin\Desktop\MicroGrid_5march\MbolaN.shp"
            #inputShapeFile=r"C:\Documents and Settings\Selin\Desktop\sampleDir\Pampaida\PampaidaPts.shp"
            #outputDir = argv[2]
            outputDir=r"C:\Documents and Settings\Selin\Desktop\MicroGrid_5march\MbolaN"
            #algorithm = argv[3]
            searchRadius = "300 meters"
        except IndexError:
            raise Error("Not enough arguments provided to script.")
        startTime = time.time()
      
        nodesByClusterID,clusterByNode,nodes,centers,LVCostDict=generateDictsFromShp(inputShapeFile,outputDir)
        
        totalCost,tree,centers,nodesByClusterID=run(centers,nodesByClusterID,clusterByNode,searchRadius)
        
        
        statsFile1= outputDir + os.sep + "ClusterSize.txt"
        statsFile2= outputDir + os.sep + "Centers.txt" 
        writeCentersToText(statsFile1,statsFile2, centers) 
        
        TreeLV=network.Network()

        TreeWithLV=addLVSeg(TreeLV,centers,nodesByClusterID)

        fileRW.genShapefile(TreeWithLV, outputDir + ".prj", outputDir + os.sep + "LV.shp")
        finishTime=time.time()       
        print "Running Time:",finishTime-startTime
    except Usage, err:
        print >>sys.stderr, err.msg

        return 2
    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

if __name__ == "__main__":
    sys.exit(main())
