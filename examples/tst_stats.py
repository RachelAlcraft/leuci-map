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
import leuci_map.mapsmanager as mapss

########## INPUTS #################
central = v3.VectorThree(1,2,3)
linear = v3.VectorThree(2,2,2)
planar = v3.VectorThree(3,2,3)
#pdb_code = "6eex" #small xray
#pdb_code = "6axz" #em


########## EXAMPLE #################
mapman = mapss.MapsManager()
mapman.set_dir(DATADIR)
pdb_codes = ["3u7z","7uly","1ejg","6eex","3nir"] 
fofcs = [(1,-1),(2,-1)]
sds = [0]
    
for pdb_code in pdb_codes:    
    t0 = datetime.datetime.now()
    mload = mapman.get_or_create(pdb_code,file=1,header=1,values=1)    
    print("------0",pdb_code,"Map:",mload.mobj.resolution,mload.mobj.exp_method,"------")
    if mload.exists():        
        for sd in sds:                        
            for fo,fc in fofcs:            
                mfunc = mfun.MapFunctions(pdb_code,mload.mobj,mload.pobj, "linear",as_sd=sd,log_level=1,fo=fo,fc=fc)
        
    print("time(s)",pdb_code,str(round((datetime.datetime.now()-t0).total_seconds(),3)))
print("~")    


            
            



                
###################################################

