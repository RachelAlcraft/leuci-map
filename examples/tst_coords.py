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
pdb_code = "1yk4"

########## EXAMPLE #################
def query_pdb():
    print("Showing pdb details", pdb_code)
    po = moad.MapLoader(pdb_code, directory=DATADIR, cif=False)
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

        central = 'A:54@FE.A'
        linear = 'A:28@CA.A'
        planar = 'A:28@O.A'
        coordsC = my_pdb.get_coords_key(central)
        coordsL = my_pdb.get_coords_key(linear)
        coordsP = my_pdb.get_coords_key(planar)
        print(coordsC,coordsL,coordsP)
        
        #print(my_pdb.lines)

        import leuci_map.mapfunctions as mfun
        import leuci_map.mapsmanager as mman

        ml = mman.MapsManager().get_or_create(pdb_code,file=1,header=1,values=1)
        mf = mfun.MapFunctions(pdb_code,ml.mobj,ml.pobj,"linear")

        cc = v3.VectorThree().from_coords(ml.pobj.get_coords_key(central))
        ll = v3.VectorThree().from_coords(ml.pobj.get_coords_key(linear))
        pp = v3.VectorThree().from_coords(ml.pobj.get_coords_key(planar))
        print(cc.A,cc.B,cc.C)
        print(ll.A,ll.B,ll.C)
        print(pp.A,pp.B,pp.C)

            
            



                
###################################################
query_pdb()
