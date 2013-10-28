#!/usr/bin/env python
"""
Contains cluster classes as well as algorithms for clustering data (kmeans / ME
from Roy).

Alex Zvoleff, aiz2101@columbia.edu
"""

import sys
import network
import kruskals

class ClusterList(list):
    def numClusts(self):
        return len(self)

    def getNodes(self):
        nodeList = []
        for clust in self:
            try:
                nodeList.extend(clust.getNodes())
            except:
                nodeList.append(clust.getNodes())
        return nodeList

    def getEdges(self):
        segList = []
        for clust in self:
            try:
                segList.extend(clust.getSegs())
            except TypeError:
                segList.append(clust.getSegs())
        return segList

    def bestInterhDist(self):
        min = None
        for clust in self:
            interDist = clust.getTotalEdgeWeight()/clust.getTotalNodeWeight()
            if interDist < min or min == None:
                min = interDist 
        return min

    def maxNodeWeight(self):
        max = None
        for clust in self:
            nodeWeight = clust.getTotalNodeWeight()
            if nodeWeight > max or max == None:
                max = nodeWeight
        return max

    def getHeaviestClust(self):
        max = None
        for clust in self:
            nodeWeight = clust.getTotalNodeWeight()
            if nodeWeight > max or max == None:
                max = nodeWeight
                maxClust = clust
        return maxClust

    def heaviestClustStats(self):
        return self.getHeaviestClust().getInterhDistance()

class MergeInfo():
    'Records indices and joining segment for two clusters to be merged.'
    def __init__(self, index1=None, index2=None, seg=None):
        if index1 > index2:
            tmp = index1
            index1 = index2
            index2 = tmp
        self._minIndex = index1
        self._maxIndex = index2
        self._seg = seg

    def __cmp__(self, other):
        if self._maxIndex > other._maxIndex:
            return 1
        elif self._maxIndex < other._maxIndex:
            return -1
        else:
            return 0

    def __str__(self):
        minIndex, maxIndex = self.getIndices()
        seg = self.getSeg()
        return "(%(minIndex)s, %(maxIndex)s, %(seg)s)" %vars()

    def __repr__(self):
        minIndex, maxIndex = self.getIndices()
        seg = self.getSeg()
        return "mergeInfo(%(minIndex)r, %(maxIndex)r, %(seg)r)" %vars()

    def getIndices(self):
        return self._minIndex, self._maxIndex

    def getSeg(self):
        return self._seg

class Cluster():
    'Stores a set of nodes, their emanating segments, and an MST connecting them.'
    def __init__(self, nodes=[], potentialSegs=set(), net=[]):
        self._nodes = nodes
        self._potentialSegs = potentialSegs
        self._net = net
        self._totalEdgeWeight = 0
        for seg in net:
            self._totalEdgeWeight += seg.getWeight()
        self._totalNodeWeight = 0
        for node in nodes:
            self._totalNodeWeight += node.getWeight()
        self._mstFlag = False

    def __str__(self):
        nodes = self.getNodes()
        return "(%(nodes)s)" %vars()

    def mergeCluster(self, other, seg):
        'Merges this cluster with another, using the given joining seg.'
        self._nodes.extend(other.getNodes())
        self._potentialSegs.update(other.getSegs())
        self._net.append(seg)
        self._net.extend(other.getNet())
        self._totalEdgeWeight += other.getTotalEdgeWeight() + seg.getWeight()
        self._totalNodeWeight += other.getTotalNodeWeight()
        self._mstFlag = False
        return self

    def splitCluster(self, seg, segDict):
        self._net.remove(seg)
        tree = network.Network()
        for seg in self._net:
            tree.addSeg(seg)
        cluster1, cluster2 = networkToClusterList(tree, segDict)
        return cluster1, cluster2

    def rebuildMST(self):
        if self.getMSTFlag() is False:
            tree = kruskals.kruskalsAlg(list(self.getSegs()),
                    self.numNodes(), None, self.getNodes())
            self.setNet(tree.getEdges())
            self.setMSTFlag()

    def getNodes(self):
        return self._nodes

    def numNodes(self):
        return len(self._nodes)

    def getSegs(self):
        return self._potentialSegs

    def numSegs(self):
        return len(self._potentialSegs)

    def getNet(self):
        return self._net

    def setNet(self, net):
        self._net = net
        self._totalEdgeWeight = 0
        for seg in net:
            self._totalEdgeWeight += seg.getWeight()

    def getTotalEdgeWeight(self):
        return self._totalEdgeWeight

    def getTotalNodeWeight(self):
        return self._totalNodeWeight

    def getInterhDistance(self):
        return self.getTotalEdgeWeight() / self.getTotalNodeWeight()

    def getMSTFlag(self):
        return self._mstFlag

    def setMSTFlag(self):
        self._mstFlag = True

    def calcMergedStats(self, other, seg):
        combinedSegWeight = self.getTotalEdgeWeight() \
                + other.getTotalEdgeWeight() + seg.getWeight()
        combinedNodeWeight = self.getTotalNodeWeight() \
                + other.getTotalNodeWeight()
        return combinedSegWeight, combinedNodeWeight

def clusterListToNetwork(clusterList):
    tree = network.Network()
    for clust in clusterList:
        if clust.numSegs() < 1:
            continue
        for seg in clust.getNet():
            tree.addSeg(seg)
    return tree 

def networkToClusterList(network, segDict):
    clustList = ClusterList()
    for netID in network.listNetIDs():
        nodes = network.getSubnetNodes(netID)
        edges = network.getSubnetEdges(netID)
        potentialSegs = set()
        for node in nodes:
            potentialSegs.update(segDict[node.getID()])
        clustList.append(Cluster(nodes, potentialSegs, edges))
    return clustList
