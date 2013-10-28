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
#import CMST_dfs
import gc
import collections
#import scipy
import csv
#import pylab
#!/usr/bin/env python
import numpy
import math
#import batchPrimsforTransformers
#from heapq import heappush, heappop
from osgeo import ogr
from modules import network
#from modules import fileRW
#from modules import prims
from collections import defaultdict

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Error(Exception):
    def __init__(self, msg):
        self.msg = msg

class Node:
    'Defines a node class, with ID, x, y, and weight attributes'
    def __init__(self, value1=None, value2=None, value3=None, value4=None, value5=None):
        self._id = value1
        self._x = value2
        self._y = value3
        self._weight = value4 # population of the settlement
        self._name=value5

    def __hash__(self):
        return hash(str(self.getID()))



    def __repr__(self):
        id = self.getID()
        x = self.getX()
        y = self.getY()
        weight = self.getWeight()
        name=self.getName()
        return "Node(%(id)r, %(x)r, %(y)r, %(weight)r), name" %vars()

    def __str__(self):
        id = self.getID()
        x = self.getX()
        y = self.getY()
        weight = self.getWeight()
        name=self.getName()
        return "(%(id)s, %(x)s, %(y)s, %(weight)s, name)" %vars()

    @property
    def __geo_interface__(self):
        'Provides an interface to allow conversion for use with Shapely.'
        return {'type': 'Point', 'coordinates': (self._x, self._y)}

    def getID(self):
        return self._id

    def setID(self,newID):
        self._id=newID

    def getX(self):
        return self._x

    def setXY(self,X,Y):
        self._x=X
        self._y=Y

    def getY(self):
        return self._y

    def getWeight(self):
        return self._weight

    def getName(self):
    	return self._name

    def setWeight(self,newWeight):
        self._weight=newWeight

    def getCoords(self):
        return self._x,self._y

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
        geomRef = feat.GetGeometryRef()
        x = geomRef.GetX()
        y = geomRef.GetY()
        FID = feat.GetFID()

        nameField= feat.GetFieldIndex("CU_NAME")
        weightField = feat.GetFieldIndex("TOTPOP")

        #nameField= feat.GetFieldIndex("Name")
        #weightField = feat.GetFieldIndex("Population")

        population = feat.GetField(weightField)
        nodeName = feat.GetField(nameField)


        nodes[FID] = Node(FID, x, y, population, nodeName) #Households
        feat = ptLayer.GetNextFeature()
    ds.Destroy()
    return nodes

def generateClusterDicts(nodes,coverDist):
    nodesByClusterID=collections.defaultdict(list)
    clustersByNodeID=collections.defaultdict(list)

    for startNode in nodes.values():
        for endNode in nodes.values():
            #dist = ((startNode.getX() - endNode.getX())**2 +
                    #(startNode.getY() - endNode.getY())**2)**(.5)
            dist=computeSphericalDistance((startNode.getX(),startNode.getY()), (endNode.getX(),endNode.getY()))
            if dist < coverDist:
                nodesByClusterID[startNode.getID()].append(endNode)
                clustersByNodeID[endNode.getID()].append(startNode.getID())

    return nodesByClusterID




def maxInClusterDist(centerNode,nodesByClusterID): #Returns maxDist within the cluster
    maxdist=0
    for node in nodesByClusterID[centerNode.getID()]: #uses the fact that centerID and ClusterID are same
        #dist=((centerNode.getX()-node.getX())**2+
            #(centerNode.getY()-node.getY())**2)**(.5)
        dist=computeSphericalDistance((centerNode.getX(),centerNode.getY()), (node.getX(),node.getY()))
        if dist>=maxdist:
            maxdist=dist
    return maxdist

def computeSphericalDistance(coordinates1, coordinates2): # coordinates in the format of (lat, long)
    """
    http://en.wikipedia.org/wiki/Great-circle_distance
    """
    # Define
    convertDegreesToRadians = lambda x: x * math.pi / 180
    # Load
    latitude1, longitude1 = map(convertDegreesToRadians, coordinates1)
    latitude2, longitude2 = map(convertDegreesToRadians, coordinates2)
    # Initialize
    longitudeDelta = longitude2 - longitude1
    earthRadiusInMeters = 6371010
    # Prepare
    y = math.sqrt(math.pow(math.cos(latitude2) * math.sin(longitudeDelta), 2) + math.pow(math.cos(latitude1) * math.sin(latitude2) - math.sin(latitude1) * math.cos(latitude2) * math.cos(longitudeDelta), 2))
    x = math.sin(latitude1) * math.sin(latitude2) + math.cos(latitude1) * math.cos(latitude2) * math.cos(longitudeDelta)
    # Return
    return earthRadiusInMeters * math.atan2(y, x) #in m



def findTheBiggestCluster(nodesByClusterID):
    maxID=0
    maxLen=0
    for nodeID in nodesByClusterID.keys():
        l=len(nodesByClusterID[nodeID])
        if l>maxLen:
            maxLen=l
            maxID=nodeID
    return maxID


'''def updateDicts(nodesByClusterID,maxID):

    copy=[]

    print "Baslangic"
    for key in nodesByClusterID.keys():
        print "key", key
        for nodes in nodesByClusterID[key]:
            print "values", nodes.getID()

    for node in nodesByClusterID[maxID]:
        copy.append(node.getID())
        del nodesByClusterID[node.getID()]
        #print node.getID(),"deleted"


    for nodeID in nodesByClusterID.keys():
        for n in nodesByClusterID[nodeID]:
            #if nodeID==506:
                #print "Copy",copy
                #for k in nodesByClusterID[506]:
                    #print "k",k.getID()
            for ID in copy:
                if nodeID==506:
                    print ID,"ID"

                if n.getID()==ID:
                    print "YEYEYEYE removed"
                    #print n.getID()
                    nodesByClusterID[nodeID].remove(n)

            #if nodeID==506:
             #   for l in nodesByClusterID[506]:
              #      print "l",l.getID()

    del copy

    print "Bitis"
    for key in nodesByClusterID.keys():
        print "key", key
        for nodes in nodesByClusterID[key]:
            print "values", nodes.getID()

    return nodesByClusterID


'''
def updateDicts(nodesByClusterID,maxID):

    copy=[]

    for node in nodesByClusterID[maxID]:
        copy.append(node.getID())
        del nodesByClusterID[node.getID()]


    for ID in copy:
        for nodeID in nodesByClusterID.keys():
            for n in nodesByClusterID[nodeID]:

                if n.getID()==ID:
                    nodesByClusterID[nodeID].remove(n)

    del copy
    return nodesByClusterID


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

            #inputShapeFile=r"C:\Users\Selin\Desktop\Shaky\png_multi_points.shp"
            inputShapeFile=r"C:\Users\Selin\Desktop\Shaky\Edwin\png_multi_points.shp"
            #outputDir = argv[2]
            outputDir=r"C:\Users\Selin\Desktop\Shaky\Edwin"
            coverDist = 1000 #meters

        except IndexError:
            raise Error("Not enough arguments provided to script.")

        startTime = time.time()

        print "Generating Dictionaries"
        nodes=generateNodeDictFromShp(inputShapeFile,outputDir)

        nodesByClusterID=generateClusterDicts(nodes,coverDist)
        #nodesByClusterID={1:[1,2,3,4],2:[2],3:[3,4],4:[3,4],5:[1,2,5,6],6:[6]}
        statsFile= outputDir + os.sep + "LTETETETET.csv"
        csvWriter = csv.writer(open(statsFile, 'wb'))
        outFile = open(statsFile,"w")

        while nodesByClusterID:
            namesInBiggestCluster=[]
            sumPopInBiggestCluster=0
            maxID=findTheBiggestCluster(nodesByClusterID)

            size=len(nodesByClusterID[maxID])

            for node in nodesByClusterID[maxID]:
            	namesInBiggestCluster.append(node.getName())
            	sumPopInBiggestCluster+=int(node.getWeight())
                namesInBiggestCluster.append('1')
            	sumPopInBiggestCluster+=int(1)

            #print "MaxID", maxID
            nodesByClusterID=updateDicts(nodesByClusterID,maxID)

            x=nodes[maxID].getX()
            y=nodes[maxID].getY()
            csvWriter.writerow([x, y, size,sumPopInBiggestCluster,namesInBiggestCluster])
        outFile.close()
        print "Total Running Time:",time.time()-startTime



    except Usage, err:
        print >>sys.stderr, err.msg

        return 2
    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1

if __name__ == "__main__":
    sys.exit(main())

