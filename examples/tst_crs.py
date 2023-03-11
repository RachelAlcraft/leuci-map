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
import leuci_xyz.vectorthree as v3
import leuci_xyz.spacetransform as sptr

from leuci_map import maploader as moad
import leuci_map.mapfunctions as mfun
import leuci_xyz.vectorthree as v3

########## INPUTS #################
pdb_code = "3u7z"
interp_method = "linear"
########## EXAMPLE #################
def query_pdb():
    print("Showing pdb details", pdb_code)
    po = moad.MapLoader(pdb_code, directory=DATADIR, cif=False)
    if not po.exists():
        po.download()
    po.load()
    my_pdb = po.pobj
    al,ac,ap = my_pdb.get_first_three()
    print(ac,al,ap)
    keyc = my_pdb.get_key(ac)
    keyl = my_pdb.get_key(al)
    keyp = my_pdb.get_key(ap)
        
    coordsc = my_pdb.get_coords_key(keyc)
    coordsl = my_pdb.get_coords_key(keyl)
    coordsp = my_pdb.get_coords_key(keyp)
    print(keyc,coordsc)
    print(keyl,coordsl)
    print(keyp,coordsp)
    cvc = v3.VectorThree().from_coords(coordsc)
    cvl = v3.VectorThree().from_coords(coordsl)
    cvp = v3.VectorThree().from_coords(coordsp)

    mf = mfun.MapFunctions(pdb_code,po.mobj,po.pobj,interp_method)

    print("Start xyz",cvc.get_key())
    c1_crs = mf.get_crs(cvc)
    print("Crs",c1_crs.get_key())
    c1_xyz = mf.get_xyz(c1_crs)
    print("return xyz",c1_xyz.get_key())

    

            
            



                
###################################################
query_pdb()
