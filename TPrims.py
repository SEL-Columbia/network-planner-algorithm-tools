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
import CMST_dfs
import gc
import collections
import scipy
#import pylab
import numpy
#import batchPrimsforTransformers
from heapq import heappush, heappop
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
    LVLenDict={}
    
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
        LVLenDict[FID]=0
        feat = ptLayer.GetNextFeature()
    ds.Destroy()
    return nodesByClusterID,clusterByNode,nodes,centers,LVLenDict

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

def primsAlg(segments, numNodes, firstNodeID, nodeDict):
    'Prim\'s Algorithm for finding a minimum spanning tree'

    
    tree = network.Network()
    segHeap = []

    # Find the shortest segment emanating from the node with the firstNodeID
    try:
        segs = nodeDict[firstNodeID]
    except KeyError:
        return tree

    leastWeight = None
    for seg in segs:
        if (seg.getWeight() < leastWeight) or (leastWeight == None):
            leastWeight = seg.getWeight()
            firstSeg = seg
    tree.addSeg(firstSeg)

    # Starter to algorithm
    # Add the segs emanating from the first two endpoints to the heap
    for endNode in [firstSeg.getNode1(), firstSeg.getNode2()]:
        addToHeap(segHeap, nodeDict[endNode.getID()])

    # Pick best from heap and repeat
    while tree.numNodes() < numNodes:
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


   
   
def run(centers,nodesByClusterID,clusterByNode,LVLenDict,sr,MV,LV,TCost,distFromT,outputDir): 
    

    maxLVLenghtInCluster=distFromT/5*6
    #maxLVLenghtInCluster=120
    minCenters=copy.deepcopy(centers)
    st=time.time()
    segments=generateSegments(minCenters,sr)
    print "genSeg", time.time()-st
    
    # To write total cost to a text file
    statFile= outputDir + os.sep + "TotalCost.txt"
    outFile = open(statFile,"w")
    
    #minTree,segments=kruskalsAlg(segments,centers) # can use either Kruskal or Prims
    m1=time.time()
    nodeDict=buildAssocDict(segments)
    
    minTree=primsAlg(segments, len(minCenters), 0, nodeDict) # 0 is the starting node of Prims algorithm
    MVLen=minTree.getTotalEdgeWeight()
    minTotalCost=MVLen*MV+len(centers)*TCost
    
    outFile.write("0 %(MVLen)f %(minTotalCost)f\n" %vars())
    minLVLenDict=copy.deepcopy(LVLenDict)
    
    minNodesByClusterID=copy.deepcopy(nodesByClusterID)
    minSeg=min(segments,key=lambda obj:obj.getWeight())
    #fixes tempCenterX , tempCenterY and maxDist of the first merge before going into while
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
        i+=1
        if i%20==0:
	   print i
	
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
        nodeDict=buildAssocDict(segments)
        #newTree,segments=kruskalsAlg(segments,centers) #returns sorted segments
        newTree=primsAlg(segments, len(centers), 0, nodeDict) # 0 is the starting node of prims.
	MVLen=newTree.getTotalEdgeWeight()
        TotalMVCost=MVLen*MV
        TotalTransformerCost=len(centers)*TCost
        del LVLenDict[mergingClusterID]
        gc.collect()
        segmentsCMST, LVLenDict[baseClusterID] =CMST_dfs.CMST(nodesByClusterID[baseClusterID],maxLVLenghtInCluster,centers[baseClusterID])
        LVLen=sum(LVLenDict.values())
        TotalLVCost=LVLen*LV
        newTotalCost=TotalMVCost+TotalTransformerCost+TotalLVCost
        
        outFile.write("%(LVLen)f %(MVLen)f %(newTotalCost)f\n" %vars()) #LVLen and MVlen are in meters, newTotalCost is in ($)whatever!
        if(newTotalCost<=minTotalCost):
            minNodesByClusterID=copy.deepcopy(nodesByClusterID)
            minTree=copy.deepcopy(newTree)
            minCenters=copy.deepcopy(centers)
            minLVLenDict=LVLenDict.copy()
            minTotalCost=newTotalCost
       
        # Calculate maxDist below for next graph and continue if it is less than 500

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
        
    outFile.close()
    return minTotalCost,minTree,minCenters,minNodesByClusterID,minLVLenDict


def addLVSeg(tree,centers,nodesByClusterID):#single points line from the root
    
    SegID=1000000

    for centerID in centers.keys():
        try:
            netID=tree.getNetID(centers[centerID])
        except:
            netID=0
            tree._nodesByNetID[0]=[]
            tree._network[netID]=[]
        for node in nodesByClusterID[centerID]:
            length=((node.getX()-centers[centerID].getX())**2+
                        (node.getY()-centers[centerID].getY())**2)**(.5)
            newSeg=network.Seg(SegID, node, centers[centerID], length)
            tree._netIDByNode[node] = netID
            tree._nodesByNetID[netID].append(node)
            tree._network[netID].append(newSeg)
   
    return tree


def writeLVDictToText(statsFile, Dict,LV):
    'Writes LVLenDict to a text file for batchPrimsForTransformers.py.'
    outFile = open(statsFile,"w")
    for key in Dict.keys():
            LVCost=Dict[key]*LV
            outFile.write("%(key)i %(LVCost)f\n" %vars())
    outFile.close()
    return 0

def writeCenterSizeToText(statsFile, Dict):
   
    outFile = open(statsFile,"w")
    for key in Dict.keys():
            size=Dict[key].getWeight()
	    x=Dict[key].getX()
	    y=Dict[key].getY()
            outFile.write("%(size)i %(x)f %(y)f \n" %vars())
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
           
            inputShapeFile=r"C:\Users\Selin\Desktop\sample\sample2.shp"
            
            #outputDir = argv[2]
            outputDir=r"C:\Users\Selin\Desktop\sample\sample2"
            #algorithm = argv[3]
            searchRadius =20000
            MV =10# Cost of MV per meter
            LV = 10 # Cost of LV per meter
            TCost=0 # Transformer Cost
            distFromT=2500
            
        except IndexError:
            raise Error("Not enough arguments provided to script.")

        startTime = time.time()
        print "Generating Dictionaries"
       
        nodesByClusterID,clusterByNode,nodes,centers,LVLenDict=generateDictsFromShp(inputShapeFile,outputDir)
       
        print "Run function starts..."
        timeBeforeRun=time.time()
        totalCost,tree,centers,nodesByClusterID,LVLenDict=run(centers,nodesByClusterID,clusterByNode,LVLenDict,searchRadius,MV,LV,TCost,distFromT,outputDir)
        print "Time for RUN:", time.time()-timeBeforeRun
        fileRW.genShapefile(tree, outputDir + ".prj", outputDir + os.sep + "MV.shp")
        
        statsFile1= outputDir + os.sep + "LVCostDict.txt" 
	statsFile2= outputDir + os.sep + "CenterSize.txt" 
        writeLVDictToText(statsFile1, LVLenDict,LV)
        writeCenterSizeToText(statsFile2, centers)
        #batchPrimsforTransformers.batchPrims(tree,centers,LVLenDict,outputDir)
       
        print "LVLength_run", sum(LVLenDict.values())
        MVLength=tree.getTotalEdgeWeight()
        
        MVCost=MVLength*MV
        numTransformer=len(centers)
        afterrun=time.clock()
        '''Use below if multipoint LV lines are needed'''
        maxLVLenghtInCluster=distFromT/5*6
        #maxLVLenghtInCluster=120
        try:
            netID=tree.getNetID(centers.values()[0])
        except:
            netID=0
            tree._nodesByNetID[0]=[]
            tree._network[netID]=[]
        start=time.time()
        for ID in centers.keys():
            nodesByNodeID={}
            segments,lvCost=CMST_dfs.CMST(nodesByClusterID[ID],maxLVLenghtInCluster,centers[ID])
           
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

        
        finish=time.time()
    
        #Tree=addLVSeg(tree,centers,nodesByClusterID) # Can be used when dedicated lines between Transformer and HHs.
        LVLength=tree.getTotalEdgeWeight()-MVLength
        print "LVLength_check", LVLength
        LVCost=LVLength*LV
        print "MVLength", MVLength
        print "Num Transformers", numTransformer
        print "Total LV Cost", LVCost
        print "Total MV Cost", MVCost
        transformerCost=numTransformer*TCost
        print "Transformer Cost",transformerCost
        print "Total Cost=" ,MVCost+LVCost+transformerCost
        
        fileRW.genShapefile(tree, outputDir + ".prj", outputDir + os.sep + "FinalGrid.shp")
        
        print "Total Running Time:",time.time()-startTime
        
        
    except Usage, err:
        print >>sys.stderr, err.msg

        return 2
    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

if __name__ == "__main__":
    sys.exit(main())
