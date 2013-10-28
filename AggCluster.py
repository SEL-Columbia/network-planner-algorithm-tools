#!/usr/bin/env python
"""
A heuristic for facility location algorithm, this version does not build a network between facilities
finds minimum possible number of facilities satisfying Dmax constraint
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


    
def mergeCluster(ClusterByNode,NodesByClusterID,Centers,segment): 
    center1,center2=segment.getNodes()
   
    centerX=(center1.getWeight()*center1.getCenterX()
        +center2.getWeight()*center2.getCenterX())/(center2.getWeight()+center1.getWeight())
    centerY=(center1.getWeight()*center1.getCenterY()
        +center2.getWeight()*center2.getCenterY())/(center2.getWeight()+center1.getWeight())

    
    weight=center2.getWeight()+center1.getWeight()
    baseClusterID=min(ClusterByNode[center1],ClusterByNode[center2])
    mergingClusterID=max(ClusterByNode[center1],ClusterByNode[center2])
   
    NodesByClusterID[baseClusterID].extend(NodesByClusterID.pop(mergingClusterID))
    
    Centers[baseClusterID].setXY(centerX,centerY)
    Centers[baseClusterID].setWeight(weight)
    
    del Centers[mergingClusterID] 
   
    for node in NodesByClusterID[baseClusterID]:
        ClusterByNode[node]=baseClusterID


        
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
    
    nodesByClusterID=collections.defaultdict(list)
    clusterByNode={}
    nodes={}
    centers={}
    LVCostDict={}
    
    feat = ptLayer.GetNextFeature()
    while feat is not None:
        nodeWeight = 1
        geomRef = feat.GetGeometryRef()
        x = geomRef.GetX()
        y = geomRef.GetY()
        FID = feat.GetFID()
        nodes[FID] = network.Node(FID, x, y, nodeWeight) #Households
        centers[FID]=network.Node(FID, x, y, nodeWeight) #Transformers (center of mass of the cluster)
        
        clusterByNode[nodes[FID]]=FID 
        nodesByClusterID[FID].append(nodes[FID])
        #LVCostDict[FID]=0
        feat = ptLayer.GetNextFeature()
    ds.Destroy()
    return nodesByClusterID,clusterByNode,nodes,centers

def generateSegments(centers,searchRadius): 
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
                segID+=1
    return segments




def maxInClusterDist(centerNode,nodesByClusterID): #Returns maxDist within the cluster
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
    


def kruskalsAlgForMID(segments,nodes):
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
    MID=tree.getTotalEdgeWeight()/numNodes
    
    return MID

def addToHeap(heap, newSegs):
    'Adds new segments to the segHeap.'
    for seg in newSegs:
        heappush(heap, seg)
    return heap

def buildAssocDict(segments):
    'Builds dictionary with nodeID key where values are segs from/to that node'
    segList = {}
    for seg in segments:
        node1, node2 = seg.getNodes()
        for nodeID in [node1.getID(), node2.getID()]:
            if segList.has_key(nodeID):
                segList[nodeID].append(seg)
            else:
                segList[nodeID] = [seg]
    return segList


   
def run(centers,nodesByClusterID,clusterByNode,sr,MV,LV,TCost,distFromT,outputDir): 
    
   
    st=time.time()
    segments=generateSegments(centers,sr)
    print "genSeg", time.time()-st
    
    
    minSeg=min(segments,key=lambda obj:obj.getWeight())
    
    if minSeg.getWeight()<=distFromT*2:
        maxDist=0 # can be anything less than 500
    else:
        maxDist=distFromT+10 # can be anything greater than 500
        print "NO CLUSTER POSSIBLE"
    
    tempCenter1,tempCenter2=minSeg.getNodes()
       
    tempCenterX=(tempCenter1.getWeight()*tempCenter1.getX()
        +tempCenter2.getWeight()*tempCenter2.getX())/(tempCenter2.getWeight()+tempCenter1.getWeight())
    tempCenterY=(tempCenter1.getWeight()*tempCenter1.getY()
        +tempCenter2.getWeight()*tempCenter2.getY())/(tempCenter2.getWeight()+tempCenter1.getWeight())
    i=0
    while(maxDist<=distFromT): 
      
	
        center1,center2=minSeg.getNodes()

	
        weight=center2.getWeight()+center1.getWeight()
        baseClusterID=min(clusterByNode[center1],clusterByNode[center2])
        mergingClusterID=max(clusterByNode[center1],clusterByNode[center2])
   
        nodesByClusterID[baseClusterID].extend(nodesByClusterID.pop(mergingClusterID))
    
        centers[baseClusterID].setXY(tempCenterX,tempCenterY)
        centers[baseClusterID].setWeight(weight)
    
        del centers[mergingClusterID] 
   
        for node in nodesByClusterID[baseClusterID]:
            clusterByNode[node]=baseClusterID

        segments=generateSegments(centers,sr) # generate segments for new graph
        
        
        
        
        try: # to check if there is a segment on the graph or there is only one cluster
            minSeg=min(segments,key=lambda obj:obj.getWeight()) 
	    maxDist,tempCenterX,tempCenterY=maxTempInClusterDist(minSeg,clusterByNode,nodesByClusterID)
            
	    if maxDist>distFromT:
		segments.sort(key=lambda obj:obj.getWeight())
	
		for seg in segments:
		    if seg.getWeight()>distFromT*2:
			break
		    else:
		        maxDist,tempCenterX,tempCenterY=maxTempInClusterDist(seg,clusterByNode,nodesByClusterID)    
                        if maxDist <=distFromT:
			    minSeg=seg            
			    break	
        except:
            break
   
    return centers,nodesByClusterID

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


def writeLVDictToText(statsFile, Dict,LV):
    'Writes LVCostDict to a text file for batchPrimsForTransformers.py.'
    outFile = open(statsFile,"w")
    for key in Dict.keys():
            
            LVCost=Dict[key]*LV
	    #print key, LVCost
            outFile.write("%(key)i %(LVCost)f\n" %vars())
    outFile.close()
    return 0

def writeCenterSizeToText(statsFile, Dict):
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
    return 0


   
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
            #inputShapeFile=r"C:\Documents and Settings\Selin\Desktop\art_data\ss.shp"
            inputShapeFile=r"C:\Users\Selin\Desktop\ruhiira\RuhiiraCombin.shp"
            #outputDir = argv[2]
            outputDir=r"C:\Users\Selin\Desktop\ruhiira\RuhiiraCombin" 
	    #outputDir=r"/home/selin/Desktop/real/random100"
            #algorithm = argv[3]
            searchRadius = 800
          
            distFromT=500
            
        except IndexError:
            raise Error("Not enough arguments provided to script.")

        startTime = time.time()
        print "Generating Dictionaries"
       
        nodesByClusterID,clusterByNode,nodes,centers=generateDictsFromShp(inputShapeFile,outputDir)
       
        print "Run function starts..."
        timeBeforeRun=time.time()
        #totalCost,tree,centers,nodesByClusterID,LVCostDict=run(centers,nodesByClusterID,clusterByNode,LVCostDict,searchRadius,MV,LV,TCost,distFromT,outputDir)
        centers,nodesByClusterID=run(centers,nodesByClusterID,clusterByNode,searchRadius,MV,LV,TCost,distFromT,outputDir)
        print "Time for RUN:", time.time()-timeBeforeRun

        statsFile= outputDir + os.sep + "Center.csv"
        csvWriter = csv.writer(open(statsFile, 'wb'))
        outFile = open(statsFile,"w")
        for center in centers.values():
            #print nodesByClusterID[center.getID()]
            centersDict=defaultdict (list)
            
            for node in nodesByClusterID[center.getID()]:
                centersDict[node.getID()]=node
            segments=generateSegments(centersDict,1000)
            MID=kruskalsAlgForMID(segments,centersDict) 
	    x=center.getX()
	    y=center.getY()
	    size=center.getWeight()
            del centersDict

	   
           
            csvWriter.writerow([x, y, size, MID])


            
        outFile.close()
	        
	
        
        #writeCenterSizeToText(statsFile2, centers)
     
        
        numTransformer=len(centers)
        afterrun=time.clock()
     
        print "Num Transformers", numTransformer        
        print "Total Running Time:",time.time()-startTime
        
        
    except Usage, err:
        print >>sys.stderr, err.msg

        return 2
    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

if __name__ == "__main__":
    sys.exit(main())

