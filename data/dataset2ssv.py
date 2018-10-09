""" Importer for T-Drive Dataset 

    takes t-drive directory and creates SSV
"""
import numpy as np;
import pandas as pd;
import os;
import sys;
from tqdm import tqdm;
from matplotlib import pyplot as plt;

def character(indir, outfile):
    print("Importing Character dataset")
    print("  In-Dir: %s" % indir)
    print("Out-File: %s" % outfile)
    fileset = sorted([os.path.join(indir,x) for x in os.listdir(indir) if x.startswith("file-")])
    print("==> Found %d files." % len(fileset))
    of = open(outfile, "w")
    of.write("x y id\n");
    for idx,f in enumerate(tqdm(fileset)):
        df = pd.read_table(f, sep=" ",header=None,  skipinitialspace=True)
        m = df.values
        m = np.cumsum(m, axis=0) # the dataset reports a smoothed derivative. integrate to get a spatial object
        m[:,-1] = np.ones(m.shape[0])*idx
        np.savetxt(of,m)
        



        

def usage(msg):
    print (msg)
    sys.exit(-1)
    

if __name__=="__main__":
    if (len(sys.argv) != 4):
        usage("Run with 3 parameters: type, dir, outfilename")
    
    cases = dict({
        "character":character
    })
    if sys.argv[1] not in cases:
        usage("Importer for %s not found." % sys.argv[0])
    cases[sys.argv[1]](sys.argv[2],sys.argv[3])
        
