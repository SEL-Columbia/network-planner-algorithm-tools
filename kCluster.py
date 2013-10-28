import os
import sys
import getopt
import math
import Pycluster
import numpy
import pylab
import collections
import itertools
from numpy import array
from scipy.cluster.vq import vq, kmeans, whiten, kmeans2
from modules import network
from modules import fileRW
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

features  = array([[ 1.9,2.3],
                   [ 1.5,2.5],
                   [ 0.8,0.6],
                   [ 0.4,1.8],
                   [ 0.1,0.1],
                   [ 0.2,1.8],
                   [ 2.0,0.5],
                   [ 0.3,1.5],
                   [ 1.0,1.0]])

def readNodesFromTxtForKmeans(filename):
    'Reads nodes  from a text file '
    try:
        fid = open(filename)
    except:
        print "ERROR: Could not load " + filename
        return 1
    line= fid.readline()
    nodes = []
    # For each line,
    #for line in fid:
        #nodeID, x, y, weight = line.split()
        #x = float(x)
        #y = float(y)
        #nodes.append((x,y))
    while "END" not in line:
        nodeID, x, y, weight = line.split()
        x = float(x)
        y = float(y)
        nodes.append((x,y))
        line = fid.readline()
    fid.close()

    #print nodes
    return nodes
    
def plotCluster(cluster, color):
    xs = [x[0] for x in cluster]
    ys = [x[1] for x in cluster]
    pylab.plot(xs, ys, 'o', color=color)
    # Plot mean
    mean = numpy.mean(cluster, axis=0)
    #meanx=numpy.mean(xs,axis=0)
    #meany=numpy.mean(ys,axis=0)
    pylab.plot([mean[0]],[mean[1]], 'x', color=color)
    #pylab.plot([meanx],[meany],'%sx' % color)
    #To clean the axis
    #axes = pylab.axes()
    #axes.set_xticks([])
    #axes.set_yticks([])
    

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help"])
        except getopt.GetoptError, msg:
            raise Usage(msg)
        try:
            
            nodesFile = argv[1]
            #nodesFile="C:\Users\Selin\Desktop\k-means\TibyOutput\TibyNodes.txt"
            
        except IndexError:
            raise Error("Not enough arguments provided to script.")

        nodes=readNodesFromTxtForKmeans(nodesFile)
        nodes = numpy.array(nodes)
        
        #results,assign=kmeans2(whiten(features),2,iter=20,thresh=0.0000000000000000001)
        #results,assignment=kmeans2(features,2,iter=100,thresh=0.0000000000000000000000000000000000000001)
        results = Pycluster.kcluster(array(nodes),nclusters=30,npass=50,method='m')
        assignments=results[0]
        #print results
        # Roy's verison of making cluster dict
        #clusterIDs = set(assignments)

        #clusterByClusterID = dict((clusterID, nodes[assignments == clusterID]) for clusterID in clusterIDs)
       #print clusterByClusterID[0]
        
        xs = [node[0] for node in array(nodes)]
        ys = [node[1] for node in array(nodes)]
        clusterByClusterID=collections.defaultdict(list)
        
        for x, y, clusterID in itertools.izip(xs,ys,assignments):
            #if clusterID not in clusterByClusterID:
                #clusterByClusterID[clusterID] = []
            clusterByClusterID[clusterID].append((x,y))
        #print clusterDictByID[3]    
        # print data

        #pylab.axis([-10000, 500000, -10000, 500000])
        pylab.figure()
        pylab.hold(True)
        colors=['r','b','g','c','m','y','k','w','#ff6c01','#00cd00']
        #colors=['r','b','g','c','m','y','k','w']
        #colors=['burlywood']
        #colors = 'rbgcmykw'
        for clusterID, color in itertools.izip(clusterByClusterID.keys(), itertools.cycle(colors)):
            #print clusterByClusterID[clusterID]
             plotCluster(clusterByClusterID[clusterID], color)
        print results[1], results[2]
        #pylab.savefig("C:\Users\Selin\Desktop\k-means\Tiby_kmedian.pdf")
        #pylab.show()
       
    except Usage, err:
        print >>sys.stderr, err.msg
        return 2

    except Error, err:
        print >>sys.stderr, "ERROR:", err.msg
        return 1 

if __name__ == "__main__":
    sys.exit(main())
