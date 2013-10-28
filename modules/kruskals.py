#!/usr/bin/env python
"""
Implements Kruskal's algoritm for finding a minimum spanning tree. Input
is a seg of seg class instances. Output is a spanning tree that is a network
class instance

Alex Zvoleff, aiz2101@columbia.edu
"""
import network

def kruskalsAlg(segments, numNodes, maxSegWeight= None):
    'Kruskal\'s algorithm for finding a minimum spanning tree'
    segments.sort(key=lambda obj:obj.getWeight())
    tree = network.Network()
    
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
        if (maxSegWeight is not None) and (segment.getWeight() > maxSegWeight):
            break
    return tree
