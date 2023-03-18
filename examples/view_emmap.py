"""
RSA 15/3/23

This loads and examines a map file and it's corresponding cif file given the pdb code
It will automatically check what kind of electron ddensity is available - xray or cryo em
"""

########## INPUTS #################
#pdb_code = "6axz"
#pdb_code = "3j9e"
pdb_code = "7a6a"


########## INPUTS #################
from pathlib import Path
DATADIR = str(Path(__file__).resolve().parent )+ "/data/"
CODEDIR = str(Path(__file__).resolve().parent.parent )+ "/src/"
import sys
sys.path.append(CODEDIR)

import leuci_map.mapobject as mobj
import leuci_map.maploader as moad


def show_pdb_map(pdb_code):
    print("Showing pdb map details", pdb_code)
    po = moad.MapLoader(pdb_code, directory=DATADIR, cif=False)
    if not po.exists():
        po.download()
    po.load()
    if po.em_loaded:
        print(po.mobj.pdb_code)
        print(po.mobj.em_code)
        print(po.mobj.pdb_link)
        print(po.mobj.ebi_link)
        print(po.mobj.em_link)
        #print(po.map_code)
        print(po.mobj.resolution)
        print(po.mobj.exp_method)
        #print(po.map_source)
        #print(po.map_link)
        #print(po.map_header)
        print(po.em_loaded)
        print(po.mobj.map_header)
        print(po.mobj.header_as_string)
        if not po.wait_for_load(log_level=1):
            po.load_values(diff=True)
        if po.values_loaded:
            print("Values...")
            print(len(po.mobj.values),"=",po.mobj.map_header["01_NC"] * po.mobj.map_header["02_NR"] * po.mobj.map_header["03_NS"])
            print(po.mobj.values[0])    
            print(po.mobj.values[len(po.mobj.values)-1])            
            print("FMS=",po.mobj.F,po.mobj.M,po.mobj.S)
            print("Diff Values...")            
            if len(po.mobj.diff_values) > 0:
                print(len(po.mobj.diff_values),"=",po.mobj.map_header["01_NC"] * po.mobj.map_header["02_NR"] * po.mobj.map_header["03_NS"])            
                print(po.mobj.diff_values[0])
                print(po.mobj.diff_values[len(po.mobj.diff_values)-1])            
            else:
                print("            ...none")
            print(po.mobj.em_link)
    #print(po._struc_dict["_pdbx_database_related.details"])
    print("~")
    

show_pdb_map(pdb_code)
