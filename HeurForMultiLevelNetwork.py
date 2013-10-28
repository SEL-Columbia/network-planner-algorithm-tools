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
#import pylab
import numpy
#import batchPrimsforTransformers
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


    
def mergeCluster(ClusterByNode,NodesByClusterID,Centers,segment): ###merge ettikten sonra eski centerlarin silindigine emin ol.
    center1,center2=segment.getNodes()
    # yeni clusterin centerlini bulur!
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
    
    del Centers[mergingClusterID] #yeni center ve cluster eklemek lazim
   
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
    segments=[]
    nodeCopy = centers.copy()
    #nodeCopy=copy.deepcopy(centers)
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



   
def run(centers,nodesByClusterID,clusterByNode,LVCostDict,sr): 
    segments=generateSegments(centers,sr)
    minTree,segments=kruskalsAlg(segments,centers) #tree1 is minTree at the beginning.
    minTotalCost=minTree.getTotalEdgeWeight()*25+len(centers)*5000

    minLVCostDict=copy.deepcopy(LVCostDict)
    minNodesByClusterID=copy.deepcopy(nodesByClusterID)
    minCenters=copy.deepcopy(centers)
    
    #fix tempCenterX , tempCenterY and maxDist of the first merge before going into while
    if segments[0].getWeight()<=1000:
        maxDist=0 # can be anything less than 500
    else:
        maxDist=501 # can be anything greater than 500
        print "NO CLUSTER POSSIBLE"

    tempCenter1,tempCenter2=segments[0].getNodes()
       
    tempCenterX=(tempCenter1.getWeight()*tempCenter1.getX()
        +tempCenter2.getWeight()*tempCenter2.getX())/(tempCenter2.getWeight()+tempCenter1.getWeight())
    tempCenterY=(tempCenter1.getWeight()*tempCenter1.getY()
        +tempCenter2.getWeight()*tempCenter2.getY())/(tempCenter2.getWeight()+tempCenter1.getWeight())
    i=0
    while(maxDist<=500):    

        
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
        TotalMVCost=newTree.getTotalEdgeWeight()*25
        TotalTransformerCost=len(centers)*5000
        del LVCostDict[mergingClusterID]
        gc.collect()
        segmentsCMST, LVCostDict[baseClusterID],householdsCMST =CMST_dfs.CMST(nodesByClusterID[baseClusterID],600,centers[baseClusterID])
        newTotalCost=TotalMVCost+TotalTransformerCost+(sum(LVCostDict.values()))*10
        if(newTotalCost<=minTotalCost):
            minNodesByClusterID=copy.deepcopy(nodesByClusterID)
            minTree=copy.deepcopy(newTree)
            minCenters=copy.deepcopy(centers)
            minLVCostDict=copy.deepcopy(LVCostDict)
            minTotalCost=newTotalCost

        # Calculate maxDist below for next graph and continue if it is less than 500

        try:
            segment=segments[0]
        except:
            break   
        maxDist,tempCenterX,tempCenterY=maxTempInClusterDist(segment,clusterByNode,nodesByClusterID)
           
    return minTotalCost,minTree,minCenters,minNodesByClusterID,minLVCostDict


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
    #print tree.getEdges()
    return tree


def writeLVDictToText(statsFile, Dict):
    'Writes LVCostDict to a text file for batchPrimsForTransformers.py.'
    outFile = open(statsFile,"w")
    for key in Dict.keys():
            LVCost=Dict[key]*10
            outFile.write("%(key)i %(LVCost)f\n" %vars())
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

        #if len(args) == 0:
            #raise Usage("No arguments specified.")
        try:
            #inputShapeFile = argv[1]
            #inputShapeFile=r"home/selin/Desktop/windows-desktop/TRANSFORMER/Tiby_Roofs/TibyQB2005Pts.shp"
            inputShapeFile=r"/home/selin/Desktop/TibyParts/T3.shp"
            #inputShapeFile=r"home/selin/Desktop/windows-desktop/TRANSFORMER/Tiby_Roofs/TibyQB2005Pts"
            #outputDir = argv[2]
            outputDir=r"/home/selin/Desktop/TibyParts/T3"
            #algorithm = argv[3]
            searchRadius = "1000 meters"
        except IndexError:
            raise Error("Not enough arguments provided to script.")
	print "T3"
        startTime = time.time()
        print "Generating Dictionaries"
       
        nodesByClusterID,clusterByNode,nodes,centers,LVCostDict=generateDictsFromShp(inputShapeFile,outputDir)
       
        print "Run function starts..."
        timeBeforeRun=time.time()
        totalCost,tree,centers,nodesByClusterID,LVCostDict=run(centers,nodesByClusterID,clusterByNode,LVCostDict,searchRadius)
        print "Time for RUN:", time.time()-timeBeforeRun
        fileRW.genShapefile(tree, outputDir + ".prj", outputDir + os.sep + "MV.shp")
        
        statsFile= outputDir + os.sep + "LVCostDict.txt" 
        writeLVDictToText(statsFile, LVCostDict)
        #batchPrimsforTransformers.batchPrims(tree,centers,LVCostDict,outputDir)
       
        print "LVCostDict=", sum(LVCostDict.values())
        MVLength=tree.getTotalEdgeWeight()
        MVCost=MVLength*25
        print "Total MV Cost", MVCost
        afterrun=time.clock()
        '''Use below if multipoint LV lines are needed'''
        try:
            netID=tree.getNetID(centers.values()[0])
        except:
            netID=0
            tree._nodesByNetID[0]=[]
            tree._network[netID]=[]
        start=time.time()
        for ID in centers.keys():
            nodesByNodeID={}
            segments,lvCost,households_CMST=CMST_dfs.CMST(nodesByClusterID[ID],600,centers[ID])
            #print "Nodes en son", households_CMST
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
        timeCMST=start-finish
        print "Time CMST", timeCMST
        #Tree=addLVSeg(tree,centers,nodesByClusterID)
        LVLength=tree.getTotalEdgeWeight()-MVLength
        LVCost=LVLength*10
        print "Total LV Cost", LVCost
        transformerCost=len(centers)*5000
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
