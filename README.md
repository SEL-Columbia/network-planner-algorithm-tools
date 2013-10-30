network-planner-algorithm-tools
===============================

About
-----

These are a series of python scripts used in the creation of __Network Planner__ ([https://github.com/modilabs/networkplanner](https://github.com/modilabs/networkplanner)), an electricity infrastructure prototyping framework developed by the [Sustainable Engineering Lab](http://modi.mech.columbia.edu/) (SEL), [Columbia University](http://columbia.edu/).

Status & Roadmap
----------------

The code in this repository does __not__ conform to SEL's [coding style guidelines](https://github.com/modilabs/StyleGuides) currently.

This repo was created for the initial purpose of preserving these scripts.

The eventual goal is to make future refinements, to convert this codebase into a reusable library.

Dependencies
------------

These scripts require that [python](http://python.org/) is installed, along with the following libraries:

* __scipy__ and __numpy__

    > Ubuntu & Debian

    > ```
sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose
```

    > For other distros and operating systems see: http://www.scipy.org/install.html

* __GDAL__ 

    > Ubuntu & Debian

    > ```
	sudo apt-get install python-gdal libgdal-dev
```

    > For other distros and operating systems: https://pypi.python.org/pypi/GDAL/

Catalog
-------

This is the current catalog of script files:

- __addOrder.py__
    > Uses prim's algorithm to build a network successively from a given point. Codes the segments according to the order they are added, and according to their percentile in the network.

    > Author: Alex Zvoleff 

- __AfricaLP.py__
    > Solves the LP model of Africa energy distribution system

    > Author: Ayse Selin Kocaman 

- __AggCluster.py__
    > A heuristic for facility location algorithm, this version does not build a network between facilities finds minimum possible number of facilities satisfying Dmax constraint

    > Author: Ayse Selin Kocaman 

- __AvgNearestFacility.py__
    > Two input shape files: Node file + facilities file Output : csv file which includes distance to the nearest facility and facility ID Calculates a statistic: Average weighted distance to a faciltiy

    > Author: Ayse Selin Kocaman 

- __AvgNearestGrid_polyline.py__
    > Calculates a statistic: Average weighted distance to a polyline. Works with the projected files only since shapely uses euclidean distance (i think). Weights of the nodes can be changed to population and can be read from the shapefile. Input shape files: Node file. Output: csv file which includes distance to the nearest road/grid (polyline).

    > Author: Ayse Selin Kocaman 

- __AvgNearestFacility_sphericalDist.py__
    > Two input shape files: Node file + facilities file Output : csv file which includes distance to the nearest facility and facility ID Calculates a statistic: Average weighted distance to a faciltiy

    > Author: Ayse Selin Kocaman 

- __batchPrims.py__
    > Implements the "Composite Prim's" algorithm to estimate how optimal network cost varies with penetration rate. Given a precalculated spanning network as a shapefile, the script tests an input set of starting points and compiles the results, returning the optimal cost for each penetration rate, and the node ID of the optimal starting point for each penetration rate.

    > Author: Alex Zvoleff 

- __batchPrimsforTransformers.py__
    > A variation of 'batchPrims.py' for a specific transformer network problem.

    > Author: Ayse Selin Kocaman 

- __calcTree.py__
    > Runs networking algorithm on a precalculated set of segments and node coordinates. Can use either Kruskal's or Prim's algorithms.

    > Author: Alex Zvoleff 

- __CMST_dfs.py__
    > A heuristic for Capaciated Minimum Spanning Tree

    > Author: Ayse Selin Kocaman 

- __distGen.py__
    > Processes a collection of shapefiles to generate distances between the points, and coordinates of the points. Search radius is a parameter of the script.

    > Author: Alex Zvoleff 

- __distStats.py__
    > Calculates interhousehold distance statistics for an input shapefile, taking account of the cumulative population using a weight parameter associated with each node.

    > Author: Alex Zvoleff 

- __DFS_BFD_Graph_Traversal.py__
    > Very simple depth first search and breath first graph traversal. This is based on the find_path written by Guido van Rossum. 
	> (source: http://code.activestate.com/recipes/576723-dfs-and-bfs-graph-traversal/)

    > Author: Bruce Wernick, http://code.activestate.com/recipes/users/4169952/

- __FacilityLoc.py__
    > A heuristic for facility location algorithm, this version does not build a network between facilities finds minimum possible number of facilities satisfying Dmax constraint

    > Author: Ayse Selin Kocaman 

- __FacilityLoc_StarConfig.py__
    > A heuristic for facility location algorithm, this version does not build a network between facilities finds minimum possible number of facilities satisfying Dmax constraint

    > Author: Ayse Selin Kocaman 

- __fullTreeCalc.py__
    > Processes a collection of shapefiles to generate distances and coordinates of points, and then runs a script to calculate an MST using the generated distances as possible edges.

    > Author: Alex Zvoleff 

- __HeurForMultiLevelNetwork.py__
    > A heuristic for facility location algorithm

    > Author: Ayse Selin Kocaman 

- __genSubnets.py__
    > Uses prim's algorithm to build a series of subnets of maximum length (in nodes) specified by an input parameter. The subnetIDs and order (within each subnet) are written to the input shapefile. The networks are calculated iteratively, with the first point chosen by choosing the point closest to a given starting point. After each subnet is formed, those nodes are removed from the set of input points, and the next subnet calculated.

    > Author: Alex Zvoleff 

- __healthcarePrims.py__
    > Uses Prim's Algorithm to build a network of a specified maximum number of nodes. The network starts at the node closest to a given starting point. The starting point is supplied as a single point in a separate shapefile.

    > Author: Alex Zvoleff 

- __kCluster.py__
    > Reads nodes from a text file, applies a k-means clustering algorithm and plots the results

    > Author: Roy Hyunjin Han

- __KruskalCluster.py__
    > Runs Kruskal, breask before last (ClusterNumber-1) segment.

    > Author: Alex Zvoleff 

- __maxDist.py__
    > Given an input point shapefile, calculates the maximum distance between any point and its nearest neighbor, and returns that distance.

    > Author: Alex Zvoleff 

- __microgridDirectLV_improved.py__
    > A heuristic for facility location algorithm

    > Author: Ayse Selin Kocaman 

- __MST.py__
    > Minimum Spanning Tree Algorithm

    > Author: Ayse Selin Kocaman 

- __PrimPenetrationRate.py__
    > Expands an existing grid with different penetartion rates 

    > Author: Alex Zvoleff 

- __SetCovering_WithNetwork.py__
    > A Set Covering Heuristic Algorithm Decides number and the locations of the facilities.

    > Author: Ayse Selin Kocaman 

- __SetCovering_withSphericalDistance_withPopandNames.py__
    > A Set Covering Heuristic Algorithm Decides number and the locations of the facilities.

    > Author: Ayse Selin Kocaman 

- __TKruskals.py__
    > A heuristic for facility location algorithm

    > Author: Ayse Selin Kocaman 

- __TPrims.py__
    > A heuristic for facility location algorithm

    > Author: Ayse Selin Kocaman 

- __treeTrim.py__
    > Trims a precalculated tree by removing the longest segment until the specified number of nodes is reached.

    > Author: Alex Zvoleff 

- __cluster.py__
    > Contains cluster classes as well as algorithms for clustering data (kmeans / ME from Roy).

    > Author: Alex Zvoleff 

- __fileRWarcgis.py__
    > Functions for reading coordinate and segment files from ArcGIS, and for converting generated networks into an ArcGIS coverage or shapefile. Also handles writing network statistics fields to a given shapefile

    > Author: Alex Zvoleff 

- __fileRWogr.py__
    > Functions for reading coordinate and segment files from ArcGIS, and for converting generated networks into an ArcGIS coverage or shapefile. Also handles writing network statistics fields to a given shapefile. Uses OGR to process shapefiles.

    > Author: Alex Zvoleff 

- __fileRW.py__
    > Chooses either fileRW-arcgis.py or fileRW-ogr.py to handle shapefile i/o. Returns an error if neither is available.

    > Author: Alex Zvoleff 

- __kruskals.py__
    > Implements Kruskal's algorithm for finding a minimum spanning tree. Input is a seg of seg class instances. Output is a spanning tree that is a network class instance

    > Author: Alex Zvoleff 

- __network.py__
    > Defines classes to allow creation of undirected networks. Separate classes are used for the full network, the constituent line segments, and their constituent nodes.

    > Author: Alex Zvoleff 

- __prims_neighborhood.py__
    > Implements Prim's algorithm for finding a minimum spanning tree. Uses a "neighborhood" to optimize the mean distance of the earlier trees (at the expense of fully spanning trees).

    > Author: Alex Zvoleff 

- __prims.py__
    > Implements Prim's algorithm for finding a minimum spanning tree.

    > Author: Alex Zvoleff 

- __readCVSthenWRITE.py__
    > Load a series of csv files and write their contents to stdout

    > Author: Ayse Selin Kocaman 

- __shapefile_editor.py__
    > Writes a field (provided as a dictionary by FID) to a shapefile.

    > Author: Alex Zvoleff 
