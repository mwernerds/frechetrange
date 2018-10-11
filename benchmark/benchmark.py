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

class cfg:
    iter=100 # 1000
    plot=False

class wrapper:
    def __init__(self, idx, t, distance_threshold):
        self.t = t
        self.distance_threshold = distance_threshold
        self.idx = idx
    def __call__(self):
        self.idx.query(self.t, self.distance_threshold)
    


if __name__=="__main__":

    distance_threshold = float(sys.argv[2]);
    
    print("Import working")

    data = pd.read_csv(sys.argv[1],sep=" ")
    query = np.random.choice(data["id"])
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
    result = idx.query(t, distance_threshold)

    w = wrapper(idx,t,distance_threshold)
    res = timeit.timeit(w, number=cfg.iter)
    print("Average query time: %f" % (res/cfg.iter))
    
    print("Query with result materialization took  %fs"%(time.time()-ts))
    print("Querying %d gives %d results" % (query, len(result)))

    if cfg.plot:
        _ = trajectories.agg(lambda x: plt.plot(x["x"],x["y"],color="gray")) #plot data (in gray)
        # Overlay results (in red)
        for q in result:
            plt.plot(q[:,0],q[:,1],color="red")
        # Draw query on top (in blue)
        plt.plot(t[:,0],t[:,1],color="blue")
        plt.show()


    ############## Second Index
    # Prepare the Index Data Structures
    idx = frechet.BBIndex();
    # Group trajectories and plot / add them individually
    _ = trajectories.agg(lambda x: idx.add(x)) #add

    # Computing the Index is not needed, it is done on insertion
    # Query
    result = idx.query(t, distance_threshold)
    w = wrapper(idx,t,distance_threshold)
    res = timeit.timeit(w, number=cfg.iter)
    print("Average query time: %f" % (res/cfg.iter))
    print("Querying %d gives %d results" % (query, len(result)))

    if cfg.plot:
        _ = trajectories.agg(lambda x: plt.plot(x["x"],x["y"],color="gray")) #plot data (in gray)
        # Overlay results (in red)
        for q in result:
            plt.plot(q[:,0],q[:,1],color="red")
        # Draw query on top (in blue)
        plt.plot(t[:,0],t[:,1],color="blue")
        plt.show()
    

    ############## Third Index
    # Prepare the Index Data Structures
    idx = frechet.TUEIndex();
    # Group trajectories and plot / add them individually
    trajectories = data.groupby(["id"]) # group
    _ = trajectories.agg(lambda x: idx.add(x)) #add

    idx.build_index();
    # Query
    result = idx.query(t, distance_threshold)
    w = wrapper(idx,t,distance_threshold)
    res = timeit.timeit(w, number=cfg.iter)
    print("Average query time: %f" % (res/cfg.iter))

    print("Querying %d gives %d results" % (query, len(result)))
    dec = frechet.FrechetDecider()
    result = [ q for q in result if dec.decide_dv(t,q, distance_threshold)]
    print("[Refined] Querying %d gives %d results" % (query, len(result)))

    if cfg.plot:
        _ = trajectories.agg(lambda x: plt.plot(x["x"],x["y"],color="gray")) #plot data (in gray)
        # Overlay results (in red)
        for q in result:
            plt.plot(q[:,0],q[:,1],color="red")
        # Draw query on top (in blue)
        plt.plot(t[:,0],t[:,1],color="blue")
        plt.show()


    
