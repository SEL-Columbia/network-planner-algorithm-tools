#!/usr/bin/env python
"""
input shape files: Node file 
Output : csv file which includes distance to the nearest road/grid (polyline)
Calculates a statistic: Average weighted distance to a polyline
Works with the projected files only since shapely uses euclidean distance (i think)
Weights of the nodes can be changed to population and can be read from the shapefile. 
Ayse Selin Kocaman
ask2170@columbia.edu
"""

import os
import sys
import time
import csv
import math
import getopt
from osgeo import ogr

from shapely.geometry import Point, LineString



class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Error(Exception):
    def __init__(self, msg):
        self.msg = msg

class Node:
    'Defines a node class, with ID, x, y, and weight attributes'
    def __init__(self, value1=None, value2=None, value3=None, value4=None):
        self._id = value1
        self._x = value2
        self._y = value3
        self._weight = value4

    def __hash__(self):
        return hash(str(self.getID()))

    def __repr__(self):
        id = self.getID()
        x = self.getX()
        y = self.getY()
        weight = self.getWeight()
        return "Node(%(id)r, %(x)r, %(y)r, %(weight)r)" %vars()

    def __str__(self):
        id = self.getID()
        x = self.getX()
        y = self.getY()
        weight = self.getWeight()
        return "(%(id)s, %(x)s, %(y)s, %(weight)s)" %vars()

    def __lt__(self, other):
        if self.getWeight() < other.getWeight():
            return True
        else:
            return False

    def __le__(self, other):
        if self.getWeight() <= other.getWeight():
            return True
        else:
            return False

    def __eq__(self, other):
        if self.getID() == other.getID():
            return True
        else:
            return False

    def __ne__(self, other):
        if self.getID() != other.getID():
            return True
        else:
            return False
    
    def __gt__(self, other):
        if self.getWeight() > other.getWeight():
            return True
        else:
            return False

    def __ge__(self, other):
        if self.getWeight() >= other.getWeight():
            return True
        else:
            return False

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
    
    nodes={}
    points={}  
    feat = ptLayer.GetNextFeature()
    while feat is not None:
        geomRef = feat.GetGeometryRef()
        X = geomRef.GetX() 
        Y = geomRef.GetY()
        '''
        Xm=feat.GetFieldIndex("x_meter")
        Ym=feat.GetFieldIndex("y_meter")
        X = feat.GetField(Xm)
        Y = feat.GetField(Ym)
        '''
        FID = feat.GetFID()
        #weightField = feat.GetFieldIndex("pop") for now only
        
        #nodeWeight = feat.GetField(weightField)
        nodeWeight=1
        nodes[FID] = Node(FID, X, Y, nodeWeight) #Creates a node dict. A node can be a HH, popCenter etc.
        points[FID]=Point([X, Y])
        feat = ptLayer.GetNextFeature()
    ds.Destroy()
    return nodes,points



def computeSphericalDistance(node1, node2):
    """
    http://en.wikipedia.org/wiki/Great-circle_distance
    """
    # Define
    convertDegreesToRadians = lambda x: x * math.pi / 180
    # Load
    longitude1, latitude1 = map(convertDegreesToRadians, node1.getCoords())
    longitude2, latitude2 = map(convertDegreesToRadians, node2.getCoords())
    # Initialize
    longitudeDelta = longitude2 - longitude1
    earthRadiusInMeters = 6371010
    # Prepare
    y = math.sqrt(math.pow(math.cos(latitude2) * math.sin(longitudeDelta), 2) + math.pow(math.cos(latitude1) * math.sin(latitude2) - math.sin(latitude1) * math.cos(latitude2) * math.cos(longitudeDelta), 2))
    x = math.sin(latitude1) * math.sin(latitude2) + math.cos(latitude1) * math.cos(latitude2) * math.cos(longitudeDelta)
    # Return
    return earthRadiusInMeters * math.atan2(y, x)



def readNetFromShpCreateLineStrings(shapefile):
    'Reads segs and nodes from the given shapefile'
    ds = ogr.Open(shapefile)
    layer = ds.GetLayer(0)
    #net = network.Network()
    feat = layer.GetNextFeature()
    nodeID=0
    lineArray=[]
    while feat is not None:
        geomRef = feat.GetGeometryRef()
        
        endPts = []
        
        for pointIndex in [0, geomRef.GetPointCount()-1]:
           
            x, y = geomRef.GetPoint(pointIndex)[0], geomRef.GetPoint(pointIndex)[1]
            
            endPts.append((x,y))
             
        lineArray.append(LineString([[endPts[0][0], endPts[0][1]], [endPts[1][0], endPts[1][1]]]))
        feat = layer.GetNextFeature()
       
    ds.Destroy()
    
    return lineArray



def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.GetoptError, msg:
            raise Usage(msg)
        try:
        
            nodesShapeFile = argv[1]
            gridFile=argv[2]
            #nodesShapeFile = "/home/selin/Desktop/trial/MbolaPts.shp"
            #nodesShapeFile = "/home/selin/Desktop/nga/nga_settlement_points.shp"
            #facilityShapeFile=argv[2]
            
            #gridFile="/home/selin/Desktop/nga/road_simplify_pts_tarred.shp"
            
            #nodesShapeFile=r"C:\Users\Selin\Desktop\Nigeria\nga_water_points_for_lga_kachia_inMeters.shp"
            #facilityShapeFile=r"C:\Users\Selin\Desktop\Nigeria\nga_health_facilities_for_lga_kachia_inMeters.shp"
            
            outputDir = argv[3]
            #outputDir="/home/selin/Desktop/nga"
                
        except IndexError:
            raise Error("Not enough arguments provided to script.")

        statsFile= outputDir + os.sep + "AvgNearestFac.csv"
        csvWriter = csv.writer(open(statsFile, 'wb'))
        outFile = open(statsFile,"w")
        #csvWriter.writerow(["Average Weighted Nearest Distance to a Facility=", avgNearestDist,"meters"])
        csvWriter.writerow(["FID","X (lat)", "Y (long)", "DistToNearestGrid",])
        
        
        
        lineArray=readNetFromShpCreateLineStrings(gridFile)
        nodes,points=generateNodeDictFromShp(nodesShapeFile,outputDir)
        
        shortestDistDict={}
        for key in points.keys():
        	shortestDistDict[key]=100000000000000000000000000000000000
        	for line in lineArray:
        		distance = line.interpolate(line.project(points[key])).distance(points[key]) 
        		if distance<=shortestDistDict[key]:
        			shortestDistDict[key]=distance
        			
        
        sumWeightedShortestDist=0
        
        for key in shortestDistDict.keys():
        	 csvWriter.writerow([key,nodes[key].getX(), nodes[key].getY(), shortestDistDict[key]])
        	 sumWeightedShortestDist=sumWeightedShortestDist+shortestDistDict[key]*nodes[key].getWeight()
        	 
        AvgWeightedShortestDist=sumWeightedShortestDist/len(nodes)
        
        print AvgWeightedShortestDist
        
        outFile.close() 
        
        
              
       
    except Usage, err:
        print >>sys.stderr, err.msg

        return 2
    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

if __name__ == "__main__":
    sys.exit(main())

