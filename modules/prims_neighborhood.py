#!/usr/bin/env python
"""
Implements Prim's algorithm for finding a minimum spanning tree. Uses a
"neighborhood" to optimize the mean distance of the earlier trees (at the
expense of fully spanning trees).

Alex Zvoleff, aiz2101@columbia.edu
"""

import sys
import heapq

import network

#numSmallest is the number of points with the smallest seg lengths to consider
#neighborhood is how far out from those segments to search
def primsAlg(segments, numNodes, firstSegID = 0, numSmallest = 0, neighborhood = 0):
    lastNumNodes = 0
    #print "Needed nodes:", numNodes
    tree = network.Network()
    nodeDict = buildAssocDict(segments)
    segHeap = []

    # Find the shortest segment emanating from the node with the firstSegID
    segs = nodeDict[firstSegID]
    leastWeight = None
    for seg in segs:
        if (seg.getWeight() < leastWeight) or (leastWeight == None):
            leastWeight = seg.getWeight()
            firstSeg = seg
    #print "FirstSeg"
    #print firstSeg

    # add the first segment to the tree
    tree.addSeg(firstSeg)
    #print tree.getEdges()

    # add the segs emanating from the first two endpoints to the heap
    for endPt in [firstSeg.getP1(), firstSeg.getP2()]:
        newSegs = nodeDict[endPt]
        addToHeap(segHeap, newSegs, tree)

    while tree.numNodes() < numNodes:
        # list the numSmallest shortest segs
        segs = heapq.nsmallest(numSmallest, segHeap)
        #print "top"
        #print "currentNumNodes", lastNumNodes 
        #print "\tTree:"
        #print "\t" + str(tree.getEdges())
        #print "\tSegs:"
        #print "\t" + str(segs)
        if len(segs) == 0:
            # Tree is finished (not all nodes contained).
            # This can happen when a small-spanning tree is
            # not input in the first place.
            break
        potentialSegs = []
        lengthRatios = []
        for seg in segs:
            endPt = getFreeEndPt(seg, tree)
            #print "\tendpt result:" + str(endPt)
            # Consider the seg if it's terminal node isn't already in the
            # cluster. network.getID will return None if a node is not
            # connected.
            if endPt == False:
                # This seg is already connected (has no free endpt)
                #print "\tRemoving", seg
                segHeap.remove(seg)
            if endPt != False:
                #print "\tin if 1"
                potentialSegs.append(seg)
                #TODO: cut primsSimple function from here
                localTree = network.Network()
                localTree.addSeg(seg)
                localSegHeap = []
                newSegs = nodeDict[endPt]
                addToHeap(localSegHeap, newSegs, localTree)
                # Now examine a network with "neighborhood" number of segs
                # created with a simplified prims algorithm from this point.
                while localTree.numNodes() < (neighborhood + 1):
                    #print "\tin while"
                    try:
                        # pop the shortest seg
                        localSeg = heapq.heappop(localSegHeap)
                    except:
                        #print "\tin while: except: break"
                        break
                    localEndPt = getFreeEndPt(localSeg, tree)
                    if localEndPt != False:
                        localTree.addSeg(localSeg)
                        localNewSegs = nodeDict[localEndPt]
                        addToHeap(localSegHeap, localNewSegs, localTree)
                #return tree
                lengthRatios.append(localTree.getTotalWeight() / 
                    localTree.numNodes())

        if len(potentialSegs) > 0:
            #print "\tin if 2"
            # Remove bestSeg from the heap, and add it and its emanating
            # segs to the heap.
            bestIndex = lengthRatios.index(min(lengthRatios))
            bestSeg = potentialSegs[bestIndex]
            #print "\tbest"
            #print "\t" + str(bestSeg)
            segHeap.remove(bestSeg)
            bestSegEndPt = getFreeEndPt(bestSeg, tree)
            tree.addSeg(bestSeg)
            newSegs = nodeDict[bestSegEndPt]
            addToHeap(segHeap, newSegs, tree)
            #print "\Segheap", segHeap
            while bestSeg in segHeap:
                #print "\tRemovedBestSegAgain"
                segHeap.remove(bestSeg)
        currentNumNodes = tree.numNodes()
        if currentNumNodes != lastNumNodes:
            print currentNumNodes
            lastNumNodes = currentNumNodes
        #print "\tNumNodes in net: " + str(tree.numNodes())
        #print "\Segheap", segHeap
#        if tree.numNodes() == 6:
#            break
    return tree

def getFreeEndPt(seg, tree):
    node1 = seg.getP1()
    node2 = seg.getP2()
    p1ID = tree.getID(node1)
    p2ID = tree.getID(node2)
    if (p1ID != p2ID):
        if p1ID == None:
            freeEndPt = node1
        else:
            freeEndPt = node2
        return freeEndPt
    else:
        # Case: neither or both ends are free
        # in this script, both ends should
        # never be free, as all segs emanate
        # from a seg already in the network
        return False

def simplePrims(segs, nodeDict, firstSeg):
    #TODO: insert code
    return tree

# TODO: fix this bottleneck
def addToHeap(heap, newSegs, tree):
    treeEdges = tree.getEdges()
    for seg in newSegs:
        if seg not in treeEdges:
            heapq.heappush(heap, seg)
    return heap

# Builds a dictionary of where key is a node
# and entries are all segs from/to that node
def buildAssocDict(segments):
    segList = {}
    for seg in segments:
        node1 = seg.getP1()
        node2 = seg.getP2()
        for node in [node1, node2]:
            if segList.has_key(node):
                segList[node].append(seg)
            else:
                segList[node] = [seg]
    return segList
