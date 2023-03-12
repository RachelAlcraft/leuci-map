"""
RSA 4/2/23

This loads and examines a map file and it's corresponding cif file given the pdb code
It will automatically check what kind of electron ddensity is available - xray or cryo em
"""

## Ensure code is importaed in path
from pathlib import Path
CODEDIR = str(Path(__file__).resolve().parent.parent )+ "/src/"
import sys
DATADIR = str(Path(__file__).resolve().parent )+ "/data/"
CODEDIR = str(Path(__file__).resolve().parent.parent )+ "/src/"
sys.path.append(CODEDIR)
from pathlib import Path
import datetime
import leuci_xyz.vectorthree as v3
import leuci_xyz.spacetransform as sptr
import leuci_map.mapobject as mobj
import leuci_map.maploader as moad
import leuci_map.mapfunctions as mfun

########## INPUTS #################
central = v3.VectorThree(1,2,3)
linear = v3.VectorThree(2,2,2)
planar = v3.VectorThree(3,2,3)
pdb_code = "6eex"
width=6
samples=100
interp_method="cubic"
new_mat = True

########## EXAMPLE #################
def single_slice():    
    print("Showing pdb map details", pdb_code)
    po = moad.MapLoader(pdb_code, directory=DATADIR, cif=False)
    if not po.exists():
        po.download()
    po.load()
    if po.em_loaded:        
        print("Loading values", pdb_code)
        po.load_values()
        if po.values_loaded:
            dt1 = datetime.datetime.now()
            mf = mfun.MapFunctions(pdb_code,po.mobj,po.pobj,interp_method)
            print("Creating slice density", pdb_code)
            vals = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=0)
            #print(vals)
            print("Creating slice radient", pdb_code)
            rads = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=1)
            #print(rads)
            print("Creating slice laplacian", pdb_code)
            laps = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=2)
            #print(laps)
            dt2 = datetime.datetime.now()
            print("Time taken=",dt2-dt1)

            
            



                
###################################################
single_slice()
