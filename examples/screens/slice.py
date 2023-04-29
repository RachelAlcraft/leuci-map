
"""
RSA 16/4/2023

This mimics the work down by the Slice screen in leuci-web views.py def slice(request):

"""
include_plot = False #for timing functions we want to exclude the plot side over which we have no control
include_timer = True
### pdb code ###
pdb_code = "1ejg"
nav = "x"
### Settings ###
settings = {}
settings["width"] = 6
settings["samples"] = 50
settings["interp"] = "bspline"
settings["navigate"] = "x"
settings["deriv"] = "three"
settings["degree"] = 3
settings["fo"] = 2
settings["fc"] = -1
settings["atomdots"] = "Y"
settings["posdots"] = "Y"
settings["navdis"] = 0.1


########### Imports ################################
#import store as stor
from pathlib import Path
# phd/leuci-async/leuci-map
CODEDIR = str(Path(__file__).resolve().parent.parent.parent)+ "/src/"
import sys
DATADIR = str(Path(__file__).resolve().parent.parent )+ "/data/"
EG_DIR = str(Path(__file__).resolve().parent.parent )+ "/OfficialExamples/"
sys.path.append(CODEDIR)

from pathlib import Path
import datetime
import plotly.graph_objs as go
from plotly.subplots import make_subplots
# Ensure libraries are imported
###############################################################################
from leuci_map import maploader as moad
from leuci_map import mapfunctions as mfun
from leuci_map import mapplotter as mpl
from leuci_map import mapsmanager as mapss
from leuci_xyz import vectorthree as v3
from leuci_xyz import spacetransform as space

mapss.MapsManager().set_dir(DATADIR)


#########################################################################


### helper functions ###
def get_slice_coords(pobj,keys = [], coords = []):
    """
    Gets the cetral; linear planar coords
    """    
    a1,a2,a3 = pobj.get_first_three()
    keyl,keyc,keyp = pobj.get_key(a1),pobj.get_key(a2),pobj.get_key(a3)
    cl,cc,cp = pobj.get_coords(a1),pobj.get_coords(a2),pobj.get_coords(a3)                    
    return [keyc,keyl,keyp],[cc,cl,cp]  
        
def upload_ed_header(pdb_code):                
    try:
        mload = mapss.MapsManager().get_or_create(pdb_code,file=1,header=1,values=2)
        return mload
    except Exception as e:
        print("ERROR:\t" + "\t" + pdb_code + ' error header uploading ' + str(e) + " " + str(datetime.datetime.now())+' hours')

def upload_ed_values(pdb_code):
    try:
        mload = mapss.MapsManager().get_or_create(pdb_code,file=1,header=1,values=1)
        return mload
    except Exception as e:
        print("ERROR:\t" + "\t" + pdb_code + ' error header uploading ' + str(e) + " " + str(datetime.datetime.now())+' hours')

def get_pdbcode_and_status(interp):        
    try:        
        mload = mapss.MapsManager().get_or_create(pdb_code,file=1,header=1,values=1)
        #mfunc = mfun.MapFunctions(pdb_code,mload.mobj,mload.pobj, "linear") #the default method is linear
        mfunc = mload.get_or_make_func(interp)                        
        return pdb_code, "x",True, True,True,mfunc
    except Exception as e:
        print("ERROR:\t" + "\t" + pdb_code + ' error header uploading ' + str(e) + " " + str(datetime.datetime.now())+' hours')


def upload_ed(pdb_code):
    try:
        upload_ed_header(pdb_code)
        upload_ed_values(pdb_code)
    except Exception as e:
        print("ERROR:\t" + "\t" + pdb_code + ' error header uploading ' + str(e) + " " + str(datetime.datetime.now())+' hours')


def runExploreScreen(interp):    
    settings["interp"] = interp
    mload = mapss.MapsManager().get_or_create(pdb_code,file=1,header=1,values=1)
    mfunc = mload.get_or_make_func(interp)                                
    mobj = mfunc.mobj
    pobj = mfunc.pobj
    if len(mobj.values) > 0:
        try:                
            [keyc,keyl,keyp],[cc,cl,cp] = get_slice_coords(pobj)                
            # we need to know some of the settings, the interp and deriv
            width = int(settings["width"])
            samples = int(settings["samples"])
            interp = settings["interp"]
            deriv = settings["deriv"]
            degree = settings["degree"]
            fo = settings["fo"]
            fc = settings["fc"]                
            central = v3.VectorThree().from_coords(cc)
            linear = v3.VectorThree().from_coords(cl)
            planar = v3.VectorThree().from_coords(cp)                                                                                                                                                                        
            dens = ["nearest"]
            rads = ["linear"]
            derivs = [0] #density
            if deriv == "three" and settings["interp"] in dens:                                    
                    derivs = [0]
            elif deriv == "three" and settings["interp"] in rads:                    
                    derivs = [0,1]
            elif deriv == "three":
                    derivs = [0,1,2]
            elif deriv == "radient":
                derivs = [1]
            elif deriv == "laplacian":
                derivs = [2]
            valss = mfunc.get_slices(central,linear,planar,width,samples,interp,derivs,fo=fo,fc=fc,log_level=1,degree=degree)
            strvals = ""
            for vls in valss:
                strvals += str(vls[0][0]) + " , "
            print("Vals returned", strvals)                                                                                                                                                                                                                                                                                                                                
        except Exception as e:                              
            print(str(e))
                            
###################################################################
def runExploreScreens():
    #tms = ["linear","linear","bspline","bspline","cubic","cubic","nearest","nearest"]
    tms = ["bspline","bspline","bspline","bspline","bspline"]
    ts = datetime.datetime.now()
    tss = []
    tse = []
    for tm in range(len(tms)):
        print(tm,tms[tm],"-------------------------")
        tss.append(datetime.datetime.now())
        runExploreScreen(tms[tm])
        tse.append(datetime.datetime.now())

    print("each time")
    for t in range(len(tms)):        
        print(tms[t],str(tse[t]-tss[t]))    
    print("Time taken to start check on file",datetime.datetime.now()-ts)


###################################################################
## PROFILER ##
import cProfile
if include_timer:
    cProfile.run('runExploreScreens()')
else:
    runExploreScreens()
    