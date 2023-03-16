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
#pdb_code = "6eex" #small xray
pdb_code = "6axz" #em
fo = 2
fc = -1
width=6
samples=100
interp_method="linear"
new_mat = True
deriv=0

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
            fo,fc = 2,-1
            print("Creating slice", pdb_code, deriv, fo, fc)
            vals = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=0,fo=fo,fc=fc)            
            fo,fc = 1,-1
            print("Creating slice", pdb_code, deriv, fo, fc)
            vals = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=1,fo=fo,fc=fc,log_level=1)
            #print(rads)
            fo,fc = 1,0
            print("Creating slice", pdb_code, deriv, fo, fc)
            vals = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=2,fo=fo,fc=fc,log_level=1)
            fo,fc = 0,1
            print("Creating slice", pdb_code, deriv, fo, fc)
            vals = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=2,fo=fo,fc=fc)
            
            dt2 = datetime.datetime.now()
            print("Time taken=",dt2-dt1)

            
            



                
###################################################
single_slice()
