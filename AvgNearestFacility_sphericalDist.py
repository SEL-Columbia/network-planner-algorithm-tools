#!/usr/bin/env python
"""
Two input shape files: Node file + facilities file
Output : csv file which includes distance to the nearest facility and facility ID
Calculates a statistic: Average weighted distance to a faciltiy
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
        weightField = feat.GetFieldIndex("pop")
        nodeWeight = feat.GetField(weightField)
        nodes[FID] = Node(FID, X, Y, nodeWeight) #Creates a node dict. A node can be a HH, popCenter etc.
        feat = ptLayer.GetNextFeature()
    ds.Destroy()
    return nodes

def generateFacDictFromShp(shapeFile,outputPath): 
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
    
    facilities={}
       
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
        weight=1
        facilities[FID] = Node(FID, X, Y, weight) #Creates a facilities dict.
        feat = ptLayer.GetNextFeature()
    ds.Destroy()
    return facilities

def generateClusterDicts(nodes, facilities): 
    nearestFacDistDict={}
    nearestFacIDDict={}
    for node in nodes.values():
        nearestDist=1000000000000000000000000000# something big
        for fac in facilities.values():
            #dist = ((node.getX() - fac.getX())**2 + 
                    #(node.getY() - fac.getY())**2)**(.5)
            dist=computeSphericalDistance(node, fac)
            if dist<=nearestDist:
                nearestDist=dist
                nearestFacID=fac.getID()
        nearestFacDistDict[node.getID()]=nearestDist
        nearestFacIDDict[node.getID()]=nearestFacID
                       
    return nearestFacDistDict, nearestFacIDDict 



def writeFieldToShp(shapefile, nearestFacDistDict):
    'Writes a field (provided as a dictionary by FID) to a shapefile.'
    #TODO: fix to allow creation of a new field in a shapefile 
    
    field="Near"
    ds = ogr.Open(shapefile)
    layer = ds.GetLayer(0)
    feat = layer.GetNextFeature()
    fieldIndex = feat.GetFieldIndex(field)
    if fieldIndex == -1:
        fieldType = getFieldType(nearestFacDistDict.values()[0])
        fieldDefn = ogr.FieldDefn(field, fieldType)
        layer.CreateField(fieldDefn)
    while feat is not None:
        FID = feat.GetFID()
        try:
            fieldValue = nearestFacDistDict[FID]
        except KeyError:
            pass
        else:
            feat.SetField(fieldIndex, fieldValue)
            feat = layer.GetNextFeature()
    #ds=None
    ds.Destroy()
    
    return 0


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





def getFieldType(fieldValue):
    'Returns OGR field type appropriate for the given value'
    if type(fieldValue) == float:
        return ogr.OFTReal
    elif type(fieldValue) == int:
        return ogr.OFTInteger
    elif type(fieldValue) == str:
        return ogr.OFTString


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
            facilityShapeFile=argv[2]
            #nodesShapeFile=r"C:\Users\Selin\Desktop\Nigeria\nga_water_points_for_lga_kachia_inMeters.shp"
            #facilityShapeFile=r"C:\Users\Selin\Desktop\Nigeria\nga_health_facilities_for_lga_kachia_inMeters.shp"
            
            outputDir = argv[3]
            #outputDir=r"C:\Users\Selin\Desktop\Nigeria"
                
        except IndexError:
            raise Error("Not enough arguments provided to script.")

        startTime = time.time()
        
        print "Generating Dictionaries"
        
        nodes=generateNodeDictFromShp(nodesShapeFile,outputDir)
        facilities=generateFacDictFromShp(facilityShapeFile,outputDir)
        nearestDict, nearestFac=generateClusterDicts(nodes, facilities)

        
        # Calculate avg.nearestDist
        sumOfWeights=0
        weightedSum=0
        
        for ID in nearestDict.keys():
            sumOfWeights=sumOfWeights+nodes[ID].getWeight()
            weightedSum=weightedSum+nodes[ID].getWeight()*nearestDict[ID]
            
        avgNearestDist=weightedSum/sumOfWeights
        
        
        statsFile= outputDir + os.sep + "AvgNearestFac.csv"
        csvWriter = csv.writer(open(statsFile, 'wb'))
        outFile = open(statsFile,"w")
        csvWriter.writerow(["Average Weighted Nearest Distance to a Facility=", avgNearestDist,"meters"])
        csvWriter.writerow(["FID","X (lat)", "Y (long)", "DistToNearestFac", "NearestFacID", "NearestFacX (lat)","NearestFacY (long)"])
        for ID in nearestDict.keys():
            FID=ID
            x=nodes[ID].getX()
            y=nodes[ID].getY()
            distToNearestFac=nearestDict[ID]
            nearestFacID=nearestFac[ID]
            nearestFacX=facilities[nearestFacID].getX()
            nearestFacY=facilities[nearestFacID].getY()
            csvWriter.writerow([FID,x, y, distToNearestFac,nearestFacID,nearestFacX,nearestFacY])
        outFile.close() 
    
          
            
        print "Average Weighted Nearest Distance to a Facility:",avgNearestDist, "meters"
        endTime=time.time()
        print"Running Time:", endTime-startTime
        
    except Usage, err:
        print >>sys.stderr, err.msg

        return 2
    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

if __name__ == "__main__":
    sys.exit(main())

