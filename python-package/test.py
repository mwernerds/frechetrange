import numpy as np;
import pandas as pd;
import frechet;
from matplotlib import pyplot as plt;
from matplotlib.path import Path
import matplotlib.patches as patches

def some_simple_and_small_tests():
    print("Create some trajectory");
    dec = frechet.FrechetDecider()
    t1 = np.array([[1,0],[2,0],[4,0],[5,0]], np.double)
    t2 = np.array([[1,0],[3,3],[5,0]])
    print(t1)
    print(t2)

    print("DV @ 2.99: %d" % (dec.decide_dv(t1,t2,2.99)))
    print("DV @ 3.01: %d" % (dec.decide_dv(t1,t2,3.01)))
    print("Decide @ 3.01: %d" % (dec.decide(t1,t2,3.01,"duetschvahrenhold")))
    print("Decide @ 3.01: %d" % (dec.decide(t1,t2,3.01,"baldusbringman")))
    # DVGrid
    grid = frechet.DVGrid();
    grid.add(t1);
    grid.add(t2);
    print(grid.as_ssv());
    grid.build_index(5)
    
    res = grid.query(t1, 2.99)
    print(res)
    res = grid.query(t1,3.01);
    print(res)



def matrix2path(m):
    print(m.values[0,])
    print(m.values.shape)
    codes = [Path.LINETO] * m.shape[0]
    codes[0] = Path.MOVETO;
    return Path(m.values,codes)
    

if __name__=="__main__":

    # Read the dataset
    sf = pd.read_csv("../data/sanfrancisco.ssv",sep=" ")
    # Prepare a query trajectory
    query = np.random.choice(sf["id"])
    t = sf.values[sf.values[:,2] == query,:]
    print("Querying for trajectory %d" % query)

    # Prepare the Index Data Structures
    grid = frechet.DVGrid();
    # Group trajectories and plot / add them individually
    trajectories = sf.groupby(["id"]) # group
    _ = trajectories.agg(lambda x: plt.plot(x["x"],x["y"],color="gray")) #plot data (in gray)
    _ = trajectories.agg(lambda x: grid.add(x)) #add

    # Compute the Index
    grid.build_index(10)
    # Query
    result = grid.query(t, 0.1)
    print("Querying %d gives %d results" % (query, len(result)))
    # Overlay results (in red)
    for q in result:
        plt.plot(q[:,0],q[:,1],color="red")
    # Draw query on top (in blue)
    plt.plot(t[:,0],t[:,1],color="blue")
    plt.show()
