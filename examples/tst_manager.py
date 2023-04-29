"""
RSA 16/4/23

Test the management of the map store class
"""

########## INPUTS #################
pdb_codes = []
pdb_codes.append("6eex")
#pdb_codes.append("6kj3")
pdb_codes.append("3j9e")
pdb_codes.append("1ejg")
#pdb_codes.append("3nir")
#pdb_codes.append("4rek")



########## INPUTS #################
from pathlib import Path
DATADIR = str(Path(__file__).resolve().parent )+ "/data/"
CODEDIR = str(Path(__file__).resolve().parent.parent )+ "/src/"
import sys
sys.path.append(CODEDIR)

import leuci_map.mapsmanager as mapss
import datetime

def show_pdb_maps():
    mapman = mapss.MapsManager()
    mapman.set_dir(DATADIR)
    runs = 3
    for r in range(runs):
        print("-----",r,"-----")
        for pdb_code in pdb_codes:
            t0 = datetime.datetime.now()
            ml = mapman.get_or_create(pdb_code,file=1,header=2,values=2)                
            print("time(s)",pdb_code,str(round((datetime.datetime.now()-t0).total_seconds(),3)))
    print("~")    
    print(mapman.print_maps())
    print("~")    


show_pdb_maps()
