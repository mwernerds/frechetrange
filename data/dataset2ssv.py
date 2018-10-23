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


def sanfrancisco(indir, outfile):
    print("Importing from a San Francisco Directory")
    print("  In-Dir: %s" % indir)
    print("Out-File: %s" % outfile)
    fileset = sorted([os.path.join(indir,x) for x in os.listdir(indir) if x.endswith(".plt")])
    print("==> Found %d files." % len(fileset))
    of = open(outfile, "w")
    of.write("x y id\n");
    for idx,f in enumerate(tqdm(fileset)):
        df = pd.read_table(f, sep=" ",header=None,  skipinitialspace=True)
        m = df.values[:,range(3)] # remove last column
        m[:,-1] = np.ones(m.shape[0])*idx
        np.savetxt(of,m)

def geolife(indir, outfile):
    print("Importing from GeoLife")
    print("  In-Dir: %s" % indir)
    print("Out-File: %s" % outfile)
#    fileset = [os.path.join(root,f) for f in files for root, _, files in os.walk(indir)]
    fileset=[]
    for root,_,files in os.walk(indir):
        fileset = fileset + [os.path.join(root, f) for f in files if f.endswith(".plt")]
    fileset = sorted(fileset)
    print("==> Found %d files." % len(fileset))
    of = open(outfile, "w")
    of.write("x y id\n");
    for idx,f in enumerate(tqdm(fileset)):
        # first skip the header
        fd = open(f, "r");
        header = [fd.readline() for _ in range(6)]
        df = pd.read_csv(fd, header=None)
        m = df.values[:,range(3)]
        m[:,-1] = np.ones(m.shape[0])*idx
        np.savetxt(of,m)
                         
#  drwxr-xr-x  3 wern_m3          1001     4096 Sep 16  2016 roma
#  drwxr-xr-x  4 wern_m3          1001     4096 Jan 10  2017 sf_large
#  drwxr-xr-x 11 wern_m3          1001     4096 Nov  5  2016 tdrive

        

def usage(msg):
    print (msg)
    sys.exit(-1)
    

if __name__=="__main__":
    if (len(sys.argv) != 4):
        usage("Run with 3 parameters: type, dir, outfilename")
    
    cases = dict({
        "character":character,
        "sanfrancisco":sanfrancisco,
        "geolife":geolife
        
    })
    if sys.argv[1] not in cases:
        usage("Importer for %s not found." % sys.argv[1])
    cases[sys.argv[1]](sys.argv[2],sys.argv[3])
        
