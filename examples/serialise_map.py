"""
RSA 4/2/23

This loads and examines a map file and it's corresponding cif file given the pdb code
It will automatically check what kind of electron ddensity is available - xray or cryo em
"""

########## INPUTS #################
#pdb_code = "6kj3"
#pdb_code = "3j9e" #em
#pdb_code = "1ejg"
pdb_code = "6eex"
#pdb_code = "3nir"
#pdb_code = "4rek"



########## INPUTS #################
from pathlib import Path
DATADIR = str(Path(__file__).resolve().parent )+ "/data/"
CODEDIR = str(Path(__file__).resolve().parent.parent )+ "/src/"
import sys
import json
import datetime
sys.path.append(CODEDIR)

import leuci_map.mapobject as mobj
import leuci_map.maploader as moad



def serialise_pdb_map(pdb_code):
    print("\n*******\nSerialising map object", pdb_code, "\n*******")
    po = moad.MapLoader(pdb_code, directory=DATADIR, cif=False)
    if not po.exists():
        po.download()
    dt1 = datetime.datetime.now()
    print("\nLoading data\t\t\t", str(dt1))
    po.load()    
    po.load_values()        
    dt2 = datetime.datetime.now()
    print("...2FoFc\t\t\t", str(dt2))
    print("\t=",str(dt2-dt1))
    po.load_values(diff=True)    
    dt2b = datetime.datetime.now()
    print("...FoFc\t\t\t\t", str(dt2b))
    print("\t=",str(dt2b-dt2))
    mob = po.mobj
    print("ONE",mob.values[0])
    dt3 = datetime.datetime.now()
    print("\nEncode into JSON formatted Data\t", str(dt3))
    #mob_json = json.dumps(mob.toJson())        
    dt4 = datetime.datetime.now()
    print("...done\t\t\t\t", str(dt4))
    print("\t=",str(dt4-dt3))
    dt5 = datetime.datetime.now()
    print("\nDecode JSON formatted Data\t", str(dt5))
    #mob2 = json.loads(mob_json)
    #mob3 = json.loads(mob2)
    #mobob = mobj.MapObject(pdb_code)
    #print(mob3["pdb_code"])
    #mobob.fromJson(mob3)
    dt6 = datetime.datetime.now()
    print("...done\t\t\t\t", str(dt6))
    print("\t=",str(dt6-dt5))
    #print(mob_json)
    #print(mobob)
    #print("TWO",mobob.npy_values[0,0,0])
    
    print("~")
     

serialise_pdb_map(pdb_code)
