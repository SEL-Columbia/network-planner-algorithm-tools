#!/usr/bin/env python
"""
Solves the LP model of Africa energy distribution system

"""

import os
import sys
import getopt
import math

import distGen
import calcTree
from modules import network
from modules import fileRW
from cvxopt import matrix
from cvxopt import solvers
from cvxopt import spmatrix
from numpy import array
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

class Node:
    'Defines a node class, with ID, x, y, and weight attributes'
    def __init__(self,value1=None,value2=None,value3=None,value4=None,value5=None,value6=None,value7=None,value8=None,value9=None,value10=None,value11=None,value12=None):
        self._id = value1
        self._x = value2
        self._y = value3
        self._demand = value4
        self._solarProd=value5
        self._windProd=value6
        self._geoProd=value7
        self._hydroProd=value8
        self._solarSupply=value9
        self._windSupply=value10
        self._geoSupply=value11
        self._hydroSupply=value12
        

    def __hash__(self):
        return hash(str(self.getID()))

    def __repr__(self):
        id = self.getID()
        x = self.getX()
        y = self.getY()
        #weight = self.getWeight()
        return "Node(%(id)r, %(x)r, %(y)r, %(weight)r)" %vars()

    def __str__(self):
        id = self.getID()
        x = self.getX()
        y = self.getY()
        #weight = self.getWeight()
        return "(%(id)s, %(x)s, %(y)s, %(weight)s)" %vars()

    @property
    def __geo_interface__(self):
        'Provides an interface to allow conversion for use with Shapely.'
        return {'type': 'Point', 'coordinates': (self._x, self._y)}

    def getID(self):
        return self._id

    def getX(self):
        return self._x

    def getDemand(self):
        return self._demand

    def getSupplyAmounts(self):
        return self._solarSupply, self._windSupply, self._geoSupply,self._hydroSupply

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
        return self._x, self._y
        
    def getSolarProdCost(self):
    	return self._solarProd
    	
    def getWindProdCost(self):
    	return self._windProd
    	
    def getGeoProdCost(self):
    	return self._geoProd
    	
    def getHydroProdCost(self):
    	return self._hydroProd
    	
    def setSolarSupply(self, newvalue):
    	self._solarSupply=newvalue
    		
    def setWindSupply(self,newvalue):
    	self._windSupply=newvalue
    	
    def setDemand(self, newvalue):
    	self._demand=newvalue
        
def readNodesFromShpforAfricaLP(shapefile):
    'Reads nodes and node weights from a point shapefile.'
    ds = ogr.Open(shapefile)
    ptLayer = ds.GetLayer(0)
    nodes = {}
    feat = ptLayer.GetNextFeature()
    while feat is not None:
        geomRef = feat.GetGeometryRef()
        FID = feat.GetFID()
        x = geomRef.GetX()
        y = geomRef.GetY()

        demandField=feat.GetFieldIndex("Dermand")
        demand=feat.GetField(demandField)

        solarProdField=feat.GetFieldIndex("Solar_Pr")
        solarProd=feat.GetField(solarProdField)
        
        windProdField=feat.GetFieldIndex("Wind_Pr")
        windProd=feat.GetField(windProdField)
        
        geoProdField=feat.GetFieldIndex("Geo_Pr")
        geoProd=feat.GetField(geoProdField)

        hydroProdField=feat.GetFieldIndex("Hydro_Pr")
        hydroProd=feat.GetField(hydroProdField)

        solarSupplyField=feat.GetFieldIndex("Solar")
        solarSupply=feat.GetField(solarSupplyField)

        windSupplyField=feat.GetFieldIndex("Wind")
        windSupply=feat.GetField(windSupplyField)

        geoSupplyField=feat.GetFieldIndex("Geo")
        geoSupply=feat.GetField(geoSupplyField)

        hydroSupplyField=feat.GetFieldIndex("Hydro")
        hydroSupply=feat.GetField(hydroSupplyField)
        
        nodes[FID] = Node(FID, x, y, demand, solarProd, windProd, geoProd, hydroProd, solarSupply, windSupply, geoSupply,hydroSupply)
        feat = ptLayer.GetNextFeature()

    
    ds.Destroy()
    return nodes

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


def generateDistanceDictionary(nodes):
    distDict = {}
    for startNode in nodes.values():
        for endNode in nodes.values():
            dist = computeSphericalDistance(startNode, endNode)
            distDict[(startNode.getID(),endNode.getID())]=dist
    return distDict



def writeSolutionToTxt(sol, outputFile):
    'Writes solution varibales to a text file.'
    ofile = open(outputFile, "w")
    
    for var in sol['x']:
    	ofile.write(" %(var).11e \n" %vars())
    
    ofile.write("END\n")
    ofile.close()
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
            #inputDir = argv[1] # Directory for dist and nodes text files.
            nodesFile="/home/selin/Desktop/AfricaLP/Africa_Updated_match21.shp"
            #nodesFile="/home/selin/Desktop/AfricaLP/3node/3nodesample.shp"
            #outputDir=argv[2] #Directory for output shapefile
            outputDir="/home/selin/Desktop/AfricaLP"
            
              
        except IndexError:
            raise Error("Not enough arguments provided to script.")


	solarTransportationCostper1000kmperGW=0.07
	windTransportationCostper1000kmperGW=0.06
	geoTransportationCostper1000kmperGW=0.03
	hydroTransportationCostper1000kmperGW=0.02
	

        numType=4
        nodes=readNodesFromShpforAfricaLP(nodesFile)
        numNodes= len(nodes)
        dist=generateDistanceDictionary(nodes)
        
        variablesFile2 = outputDir + os.sep + "variables_deneme.txt"
        
        '''
        nodes[0].setSolarSupply(1)
      	nodes[0].setWindSupply(9)
      	nodes[0].setDemand(7)
      	
      	nodes[1].setSolarSupply(3)
      	nodes[1].setWindSupply(2)
      	nodes[1].setDemand(3)
      	
      	nodes[2].setSolarSupply(1)
      	nodes[2].setWindSupply(4)
      	nodes[2].setDemand(2)
      	'''
        
        b=[]
        for node in nodes.values():
            b.append(-1*node.getDemand())# SHOULD BE MULTIPLIED BY -1
            
        
        for node in nodes.values():
            solar, wind,geo,hydro=node.getSupplyAmounts()
            
            b.append(solar)
            b.append(wind)
            b.append(geo)
            b.append(hydro)


        for i in range(0,numNodes*numType):
            b.append(0)

        for i in range(0,(numNodes*numNodes-numNodes)*numType+numType*numNodes):
            b.append(0)
        
        print len(b)
  

       

       
        
        #numNodes=5
        consts={}
        #First n constraints
        for ID in range(0,numNodes):
            consts[ID]=[]
            for i in range(0,numNodes):
                for j in range(0, numNodes):
                    if i==j:
                        continue
                    for k in range(1,numType+1):
                       
                        if i==ID:
                            consts[ID].append(1)
                            #print "ONE"
                        if j==ID:
                            consts[ID].append(-1)
                            #print "MINUS"
                        if i!=ID and j!=ID:
                            consts[ID].append(0)
                            #  print "ZERO"

            for l in range(0,numNodes):
            	for k in range(1,numType+1):
                    if l==ID:
                	consts[ID].append(-1)
                    else:
                        consts[ID].append(0)

       
        #print "LEN", len(consts[4])

        
        #4n Constraints for all xil<=Sil
        count=0
        for ID in range(numNodes,(numType+1)*numNodes): # su an 2. 4 olacak aslinda

            consts[ID]=[]
            for i in range(0,numNodes):
                for j in range(0, numNodes):
                    if i==j:
                        continue
                    for k in range(1,numType+1):
                        consts[ID].append(0)
                         
            for l in range(0,numNodes*numType):
                 if l==count:
                     consts[ID].append(1)
                 else:
                     consts[ID].append(0)

            count+=1
           
       
        #print "A"

        
        #4n Constraints for all Sumj(Yijk)<=Xik      
        count=0
        for ID in range((numType+1)*numNodes,2*numType*numNodes+numNodes): 
            consts[ID]=[]
            for i in range(0,numNodes):
               
                for j in range(0,numNodes):
                    if i==j:
                        continue
                    for k in range(1,numType+1):
                        
                        if count/numType==i and count%numType==(k-1):
                            consts[ID].append(1)
                        else:
                            consts[ID].append(0)    
                         
            for l in range(0,numNodes*numType):
                if l==count:
                    consts[ID].append(-1)
                else:
                    consts[ID].append(0)

            count+=1
            
          
            
	
        
	# variables should be positive constraints
	
        I= spmatrix(-1, range((numNodes*numNodes-numNodes)*numType+numType*numNodes), range((numNodes*numNodes-numNodes)*numType+numType*numNodes))
        
        AList=[]   

        for i in range(0,numNodes*(2*numType+1)): # *2 vardi ilk halinde
            if len(consts[i])!=(numNodes*numNodes-numNodes)*numType+numType*numNodes:
                print "ERROR in MATRIX CONSTRUCTION"

            AList.append(consts[i])
        print len(AList)
        #ATrans=matrix(AList,((numNodes*numNodes-numNodes)*numType+numType*numNodes,(2*numType+1)*numNodes),'d')
        ATrans=matrix(AList,((numNodes*numNodes-numNodes)*numType+numType*numNodes,(2*numType+1)*numNodes),'d')
        Ahalf=ATrans.trans()
  
        A=matrix([Ahalf,I])
        print A.size, "SIZE A"
            
        B=matrix(b,(((2*numType+1)*numNodes+(numNodes*numNodes-numNodes)*numType+numType*numNodes),1),'d')

        #B=matrix(b,((numType+1)*numNodes+(numNodes*numNodes-numNodes)*numType+numType*numNodes,1),'d')
        # Built cost array (c)
        c=[]
        
        for i in range(0,numNodes):
               for j in range(0, numNodes):
                   if i==j:
                       continue
                   for k in range(1,numType+1):
                       if k==1:
                           c.append(dist[(nodes[i].getID(),nodes[j].getID())]*solarTransportationCostper1000kmperGW/1000000)
                       if k==2:
                           c.append(dist[(nodes[i].getID(),nodes[j].getID())]*windTransportationCostper1000kmperGW/1000000)
                       if k==3:
                           c.append(dist[(nodes[i].getID(),nodes[j].getID())]*geoTransportationCostper1000kmperGW/1000000)
                       if k==4:
                           c.append(dist[(nodes[i].getID(),nodes[j].getID())]*hydroTransportationCostper1000kmperGW/1000000)
                                                   
                         
        for i in range(0,numNodes):
             	for k in range(1,numType+1):
       		    if k==1:
                        c.append(nodes[i].getSolarProdCost())
                    if k==2:
                        c.append(nodes[i].getWindProdCost())
                    if k==3:
                        c.append(nodes[i].getGeoProdCost())
                    if k==4:
                        c.append(nodes[i].getHydroProdCost())
         
        
      
             
        print B.size       
        C=array(c)
        C2=matrix(C,((numNodes*numNodes-numNodes)*numType+numType*numNodes, 1),'d')
	
        sol=solvers.lp(C2,A,B)
        
        variablesFile = outputDir + os.sep + "variables_march21.txt"
        indexFile=outputDir+os.sep+"index_march21.txt"
    	writeSolutionToTxt(sol, variablesFile)
        
        
        print sol['status']
        print sol['x']
        
        ofile = open(indexFile, "w")
    
    	for i in range(0,numNodes):
    		
    		for k in range(1,numType+1):
    	        	ofile.write(" %(i)i %(k)i\n" %vars())
    
        for i in c:
        	ofile.write(" %(i)f \n" %vars())
        ofile.write("END\n")
        ofile.close()
        
        
        
        
        
        
   	sum=0
   	
        for i in range(0,len(c)):
        	for l in range(0,len(sol['x'])):
        		if l==i:
        	 		
        			sum=sum+c[i]*sol['x'][l]
        		
        print "OBJECTIVE VALUE", sum
        
       

    except Usage, err:
            
        print >>sys.stderr, err.msg
        return 2

    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

if __name__ == "__main__":
    sys.exit(main())


