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

from leuci_map import pdbobject as pobj
from leuci_map import maploader as moad

########## INPUTS #################
pdb_code = "6eex"

########## EXAMPLE #################
def query_pdb():
    print("Showing pdb details", pdb_code)
    po = moad.MapLoader(pdb_code, directory=DATADIR, cif=True)
    if not po.exists():
        po.download()
    if not po.success():
        print("still not there")
    else:
        po.load()
        my_pdb = po.pobj
        a1,a2,a3 = my_pdb.get_first_three()
        print(a1,a2,a3)

        key1 = my_pdb.get_key(a1)
        key2 = my_pdb.get_next_key(key1,0)
        print("GONE BACK",key2)
        coords1 = my_pdb.get_coords_key(key1)
        coords2 = my_pdb.get_coords_key(key2)
        print(key1, key2)
        print(coords1, coords2)

        key1 = my_pdb.get_key(a2)
        key2 = my_pdb.get_next_key(key1,0)
        coords1 = my_pdb.get_coords_key(key1)
        coords2 = my_pdb.get_coords_key(key2)
        print(key1, key2)
        print(coords1, coords2)

        key1 = my_pdb.get_key(a3)
        key2 = my_pdb.get_next_key(key1)
        coords1 = my_pdb.get_coords_key(key1)
        coords2 = my_pdb.get_coords_key(key2)
        print(key1, key2)
        print(coords1, coords2)

        key_next = key1
        for i in range(20):
            key_next = my_pdb.get_next_key(key_next)
            print(i,key_next)

        key = "A:35@C.A"
        atc = my_pdb.get_atm_key(key)
        print(atc)

        key_next = my_pdb.get_next_key(key,-1)
        key_next = my_pdb.get_next_key("",-1)
        print("KEY=",key_next)    
        atc = my_pdb.get_atm_key(key_next)
        print("DIC=",atc)
        #aap = atc["aa"]              

        key = "A:35@x.A"
        atc = my_pdb.get_atm_key(key)
        print(atc)
        cc = my_pdb.get_coords_key(key)
        print(cc)
        atc = my_pdb.get_atm_key(key)                
        

            
            



                
###################################################
query_pdb()
