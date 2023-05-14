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
central = v3.VectorThree(-6.172,-15.897,3.134)
linear = v3.VectorThree(-6.447,-15.681,1.689)
planar = v3.VectorThree(-5.583,-16.994,3.507)
pdb_code = "7uly" #em
width=6
samples=10
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
            mf = mfun.MapFunctions(pdb_code,po.mobj,po.pobj,"linear")
            print("Creating neighbours", pdb_code)
            naybs = mf.get_slice_neighbours(central,linear,planar,width,samples,[0,0.5],log_level=1)
            print(naybs)            
            dt2 = datetime.datetime.now()
            print("Time taken=",dt2-dt1)

            
            



                
###################################################
single_slice()
