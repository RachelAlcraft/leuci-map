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
pdb_code = "6eex" #small xray
#pdb_code = "6axz" #em
########## EXAMPLE #################
def single_proj():    
    print("Showing pdb map details", pdb_code)
    mload = moad.MapLoader(pdb_code, directory=DATADIR, cif=False)
    if not mload.exists():
        mload.download()
    mload.load()
    if mload.em_loaded:        
        print("Loading values", pdb_code)
        mload.load_values()
        if mload.values_loaded:
            dt1 = datetime.datetime.now()
            mf = mfun.MapFunctions(pdb_code,mload.mobj,mload.pobj,"nearest")
            print("Creating slice project", pdb_code)
            valsxy = mf.get_map_projection("xy")
            valsyz = mf.get_map_projection("yx")
            valszx = mf.get_map_projection("zx")
            print(valsxy)
            print(valsyz)
            print(valszx)
            dt2 = datetime.datetime.now()
            print("Time taken=",dt2-dt1)
            # Now find the extremes of the atoms and the matching slice
            #xx,yy,zz = mf.get_atoms_extremes()
            xs,ys,zs,vs,xx,yy,zz = mf.get_atoms_projection("linear",log_level=1)
            valsxy_a = mf.get_map_projection("xy",xx[0],xx[1],yy[0],yy[1])
            valsyz_a = mf.get_map_projection("yz",yy[0],yy[1],zz[0],zz[1])
            valszx_a = mf.get_map_projection("zx",zz[0],zz[1],xx[0],xx[1])
            print(valsxy_a)
            print(valsyz_a)
            print(valszx_a)


            
            



                
###################################################
single_proj()
