"""
benchmark.py 


This file is used to (weakly) benchmark the implementations on some SSV.
It should not be used to draw conclusions about the algorithm, as some simplifications needed to 
restructure the code imply additional overheads like memory copying adding bookkeeping complexity to 
the involved implementations. If you are interested in real results regarding the performance,
rely on the C++ implementations and read them carefully in order to avoid to destroy performance
aspects of these highly optimized implementations.



"""
import sys;
import frechet;
import time;
import pandas as pd;
import numpy as np;
from matplotlib import pyplot as plt;
import sys;
import timeit;


if __name__=="__main__":
    filename = "../data/geolive.ssv"
    distance_threshold = 0.1;
    query = 5339
    result=dict();
    
    print("Import working")
    data = pd.read_csv(filename,sep=" ")
    t = data.values[data.values[:,2] == query,:]
    print("Querying for trajectory %d" % query)
    print("Phase 1: Prepare the index structure")
    trajectories = data.groupby(["id"]) # group

    # Prepare the Index Data Structures
    idx = frechet.DVGrid();
    # Group trajectories and plot / add them individually
    ts = time.time()
    _ = trajectories.agg(lambda x: idx.add(x)) #add
    print( "Copying to C++ took %fs"% (time.time()-ts));
    ts = time.time()
    # Compute the Index
    idx.build_index(distance_threshold*10)
    print("Building the index took %fs" %(time.time()-ts))
    
    # Query
    ts = time.time()
    result["dv"] = idx.query(t, distance_threshold)
    
    print("Query with result materialization took  %fs"%(time.time()-ts))

    ############## Second Index
    # Prepare the Index Data Structures
    idx = frechet.BBIndex();
    # Group trajectories and plot / add them individually
    _ = trajectories.agg(lambda x: idx.add(x)) #add

    # Computing the Index is not needed, it is done on insertion
    # Query
    result["bb"] = idx.query(t, distance_threshold)

    ############## Third Index
    # Prepare the Index Data Structures
    idx = frechet.TUEIndex();
    # Group trajectories and plot / add them individually
    trajectories = data.groupby(["id"]) # group
    _ = trajectories.agg(lambda x: idx.add(x)) #add

    idx.build_index();
    # Query
    result["tue"] = idx.query(t, distance_threshold)

    print("Refinement")
    dec = frechet.FrechetDecider()
    result["tue_refined"] = [ q for q in result["tue"] if dec.decide_dv(t,q, distance_threshold)]

    def frechet_distance(a,b, start=0.5, eps=1E-6):
#        print("Frechet Distance invocation")
        dec=frechet.FrechetDecider()
        high = start;
        low = 0
        # phase 1: grow exponentionally until decider 
        while not dec.decide_dv(a,b,high):
            print(dec.decide_dv(a,b,high))
            low = high
            high = high* 2
#            print("Grow to %f->%f" % (low, high))
        # phase 2: intervallschachtelung (lovely german word ;-)
        low = 0
        while abs(low - high) > eps:
            m = (low + high) / 2
            print(m)
            if not dec.decide_dv(a,b,m):
                low = m
            else:
                high=m
#            print("Refined to %f->%f"% (low,high))
        return (low + high) / 2

    

    tue_d = [frechet_distance(t,q) for q in result["tue"]]
    #for q in result["tue"]:
    #    print(frechet_distance(t,q))
    print(tue_d)
        
    
    

    plt.subplot(221)
    for q in result["dv"]:
        plt.plot(q[:,0],q[:,1],color="red")
    plt.plot(t[:,0],t[:,1],color="black") 
    plt.subplot(222)
    for q in result["bb"]:
        plt.plot(q[:,0],q[:,1],color="red")
    plt.plot(t[:,0],t[:,1],color="black") 
    plt.subplot(223)
    for q in result["tue"]:
        plt.plot(q[:,0],q[:,1],color="red")
    plt.plot(t[:,0],t[:,1],color="black") 
    plt.subplot(224)
    for q in result["tue_refined"]:
        plt.plot(q[:,0],q[:,1],color="red")
    plt.plot(t[:,0],t[:,1],color="black") 
    plt.show()

    
