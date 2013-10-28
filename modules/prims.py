#!/usr/bin/env python
"""
Implements Prim's algorithm for finding a minimum spanning tree.

Alex Zvoleff, aiz2101@columbia.edu
"""

import sys
from heapq import heappush, heappop

import network

def primsAlg(segments, numNodes, firstNodeID, excludedNodes = [], 
        nodeDict = None):
    'Prim\'s Algorithm for finding a minimum spanning tree'

    # How to start from old tree
        # make parameter with existing tree (instance of network.Network)
        # can use old buildAssocDict (even though it's not optimal)
        # when getting the best first segment (algorithm starter segment),
            # add a check to make sure the first segment is not in the old tree
    
    tree = network.Network()
    if nodeDict == None:
        nodeDict = buildAssocDict(segments, excludedNodes)
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

def buildAssocDict(segments, excludedNodes = []):
    'Builds dictionary with nodeID key where values are segs from/to that node'
    segList = {}
    for seg in segments:
        node1, node2 = seg.getNodes()
        if (node1 in excludedNodes) or (node2 in excludedNodes):
            continue
        else:
            for nodeID in [node1.getID(), node2.getID()]:
                if segList.has_key(nodeID):
                    segList[nodeID].append(seg)
                else:
                    segList[nodeID] = [seg]
    return segList
