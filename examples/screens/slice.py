
"""
RSA 16/4/2023

This mimics the work down by the Slice screen in leuci-web views.py def slice(request):

"""
include_plot = False #for timing functions we want to exclude the plot side over which we have no control
include_timer = False
### pdb code ###
pdb_code = "6eex"
nav = "x"

### Settings ###
settings = {}
settings["width"] = 6
settings["samples"] = 50
settings["interp"] = "bspline"
settings["central"] = ""
settings["linear"] = ""
settings["planar"] = ""
settings["keyc"] = "A:707@C.A"
settings["keyl"] = "A:707@CA.A"
settings["keyp"] = "A:707@O.A"
settings["navigate"] = "x"
settings["deriv"] = "three"
settings["degree"] = 3
settings["fo"] = 2
settings["fc"] = -1
settings["atomdots"] = "Y"
settings["posdots"] = "Y"
settings["navdis"] = 0.1


########### Imports ################################
import store as stor
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
from leuci_xyz import vectorthree as v3
from leuci_xyz import spacetransform as space


#########################################################################


### helper functions ###
def get_slice_settings(keys = [], coords = []):
    """
    Returns pdb info from current state or downloads from the ebi
    """    
    central, linear, planar = "(2.884,8.478,4.586)","(3.475,7.761,5.794)","(1.791,9.045,4.633)"
    settings["central"] = central
    settings["linear"] = linear
    settings["planar"] = planar
    return settings

def upload_ed_header(pdb_code):
    try:
        stre = stor.Store()
        if stre.exists_interper(pdb_code):
            mload,dt = stre.get_interper(pdb_code)
            return mload
        #############################################
        print("uploading header...")        
        import json        
        mload = moad.MapLoader(pdb_code, directory=DATADIR)        
        mload.load()                                                
        stre.add_loader(pdb_code,mload)
        print("added",pdb_code,"to store")                    
        return mload
    except Exception as e:
        print("ERROR:\t" + "\t" + pdb_code + ' error header uploading ' + str(e) + " " + str(datetime.datetime.now())+' hours')

def upload_ed_values(pdb_code,mload):    
    print("uploading values...")    
    try:                        
        if not mload.wait_for_load(log_level=1):
            mload.load_values(diff=True)
            stre = stor.Store()
            mfunc = mfun.MapFunctions(pdb_code,mload.mobj,mload.pobj, "linear") #the default method is linear
            stre.add_interper(pdb_code,mfunc)        
            print(mload.mobj.map_header["01_NC"],mload.mobj.map_header["02_NR"],mload.mobj.map_header["03_NS"])
            print("FMS=",mload.mobj.F,mload.mobj.M,mload.mobj.S)                    
            print("INFO:\t" + "\t" + pdb_code + ' valued loaded at '+str(datetime.datetime.now())+' hours')
    except Exception as e:
        print("ERROR:\t" + "\t" + pdb_code + ' error loading values ' + str(e) + " " + str(datetime.datetime.now())+' hours')

def get_pdbcode_and_status(ret="DATA"):
    on_file, in_loader, in_interp,mobj, mfunc = True,False,False,None,None            
    stre = stor.Store()
    in_loader = stre.exists_loader(pdb_code)
    if in_loader:
        mload,dt = stre.get_loader(pdb_code)
        mobj = mload.mobj    
    in_interp = stre.exists_interper(pdb_code)
    if in_interp:
        mfunc,dt = stre.get_interper(pdb_code)
        mobj = mfunc.mobj
    if not in_loader:
        mload = moad.MapLoader(pdb_code, directory=DATADIR, cif=False)
        on_file = mload.exists()
        if not on_file:
            mload.download()
        if on_file:
            mload = moad.MapLoader(pdb_code, directory=DATADIR)            
            mload.load()
            mobj = mload.mobj   
            on_file = mload.exists()    
    if ret == "FUNC":
        return pdb_code, nav,on_file, in_loader, in_interp,mfunc
    else:
        return pdb_code, nav,on_file, in_loader, in_interp,mobj

def upload_ed(pdb_code):
    mload = upload_ed_header(pdb_code)
    upload_ed_values(pdb_code,mload)


def runExploreScreen():
    context = {}
    print("Explore screen")    
    ts = datetime.datetime.now()    
    t1 = datetime.datetime.now()    
    pdb_code, nav,on_file, in_loader, in_interp, mfunc = get_pdbcode_and_status(ret="FUNC")
    print("page", pdb_code, nav,on_file, in_loader, in_interp)            
    make_slice = True
    print("Time taken to start check on file",datetime.datetime.now()-t1)
    t1 = datetime.datetime.now()
    if not on_file:        
        print("error not on file")
        return 0
    # upload into memory
    upload_ed(pdb_code)
    pdb_code, nav,on_file, in_loader, in_interp, mfunc = get_pdbcode_and_status(ret="FUNC")
    print("Time taken to start upload",datetime.datetime.now()-t1)
    t1 = datetime.datetime.now()                            
    if make_slice: # then it is in store and we are going to return a slice view
        mobj = mfunc.mobj
        pobj = mfunc.pobj
        if len(mobj.values) > 0:
            try:
                print("getting slice... ")
                a1,a2,a3 = pobj.get_first_three()
                keyl,keyc,keyp = pobj.get_key(a1),pobj.get_key(a2),pobj.get_key(a3)
                cl,cc,cp = pobj.get_coords(a1),pobj.get_coords(a2),pobj.get_coords(a3)                    
                # defaults from pdb if needed
                settings_dic = get_slice_settings([keyc,keyl,keyp],[cc,cl,cp])
                print("Time taken to get slice settings",datetime.datetime.now()-t1)
                t1 = datetime.datetime.now()

                # we need to know some of the settings, the interp and deriv
                width = int(settings_dic["width"])
                samples = int(settings_dic["samples"])
                interp = settings_dic["interp"]
                deriv = settings_dic["deriv"]
                degree = settings_dic["degree"]
                fo = settings_dic["fo"]
                fc = settings_dic["fc"]

                keyc = settings_dic["keyc"]                
                keyl = settings_dic["keyl"]
                keyp = settings_dic["keyp"]
                central = v3.VectorThree().from_coords(settings_dic["central"])
                linear = v3.VectorThree().from_coords(settings_dic["linear"])
                planar = v3.VectorThree().from_coords(settings_dic["planar"])                
                                                                
                nav = settings_dic["navigate"]                    
                # the key to finding the slice is the central-linear-planar coordinates                    
                if nav == "x": #then we use the given coordinates            
                    pass       
                elif nav[:2] == "A:":#if nav == "a:" then we use the givern atoms, 0 the given ones, -1 and +1 obvious                                                                                 
                    offset = int(nav[2:])                        
                    print("offset",offset)
                    if offset != 0:
                        keyc = pobj.get_next_key(keyc,offset)
                        keyl = pobj.get_next_key(keyl,offset)
                        keyp = pobj.get_next_key(keyp,offset)
                    print(keyc,keyl,keyp)
                    cc = pobj.get_coords_key(keyc)
                    cl = pobj.get_coords_key(keyl)
                    cp = pobj.get_coords_key(keyp)
                    print(cc,cl,cp)
                    central = v3.VectorThree().from_coords(cc)
                    linear = v3.VectorThree().from_coords(cl)
                    planar = v3.VectorThree().from_coords(cp)
                    print(central.get_key(), linear.get_key(), planar.get_key())
                elif nav[:2] == "N:":
                    navi = nav[2:]
                    navdis = settings_dic["navdis"]
                    print("Navigatin",navi)
                    spc = space.SpaceTransform(central, linear, planar)
                    central = spc.navigate(central,navi,navdis)
                    linear = spc.navigate(linear,navi,navdis)
                    planar = spc.navigate(planar,navi,navdis)
                                                        
                atc = pobj.get_atm_key(keyc)                
                atl = pobj.get_atm_key(keyl)                
                atp = pobj.get_atm_key(keyp)                                                
                aac,aal,aap = atc["aa"],atl["aa"],atp["aa"]  
                context['nav'] = "x"
                context["aac"] = aac
                context["aal"] = aal
                context["aap"] = aap 
                for name, val in settings_dic.items():
                    context[name] = val                                                
                context["value_check"] = mobj.values[0]
                context["value_len"] = len(mobj.values)
                context["central"] = central.get_key() 
                context["linear"] = linear.get_key()
                context["planar"] = planar.get_key()
                context["keyc"] = keyc
                context["keyl"] = keyl
                context["keyp"] = keyp
                
                # find the distance from given atoms to given coords
                cc = v3.VectorThree().from_coords(pobj.get_coords_key(keyc))
                ll = v3.VectorThree().from_coords(pobj.get_coords_key(keyl))
                pp = v3.VectorThree().from_coords(pobj.get_coords_key(keyp))
                context["disc"] = round(cc.distance(central),4)
                context["disl"] = round(ll.distance(linear),4)
                context["disp"] = round(pp.distance(planar),4)
                
                vals,rads,laps = [[]],[[]],[[]]
                context["density_mat"] = vals
                context["radient_mat"] = rads
                context["laplacian_mat"] = laps

                context["den_blocknone"] = "" # "display:block"
                context["rad_blocknone"] = ""
                context["lap_blocknone"] = ""
                context["three_blocknone"] = ""
                context["one_blocknone"] = "display:none;visibility: collapse"
                context["other_blocknone"] = "display:none;visibility: collapse"

                print("Time taken to get prepare to get slice",datetime.datetime.now()-t1)
                t1 = datetime.datetime.now()

                derivs = [0] #density
                if deriv == "three":
                    if settings_dic["interp"] == "nearest":
                        derivs = [0]
                    elif settings_dic["interp"] == "linear":
                        derivs = [0,1]
                    else:
                        derivs = [0,1,2]
                elif deriv == "radient":
                    derivs = [1]
                elif deriv == "laplacian":
                    derivs = [2]
                valss = mfunc.get_slices(central,linear,planar,width,samples,interp,derivs,fo=fo,fc=fc,log_level=1,degree=degree)

                for vls in valss:
                    print("Val returned", vls[0][0])
                                                                    
                if deriv == "three":
                    vals = valss[0].tolist()
                    if settings_dic["interp"] != "nearest":
                        rads = valss[1].tolist()
                        if settings_dic["interp"] != "linear":
                            laps = valss[2].tolist()
                    context["density_mat"] = vals
                    context["radient_mat"] = rads
                    context["laplacian_mat"] = laps
                    print(vals[0][0])
                    
                else:                        
                    context["den_blocknone"] = "display:none;visibility: collapse"
                    context["rad_blocknone"] = "display:none;visibility: collapse"
                    context["lap_blocknone"] = "display:none;visibility: collapse"
                    context["three_blocknone"] = "display:none;visibility: collapse"
                    context["other_blocknone"] = "display:none;visibility: collapse"                        
                    if deriv == "density":                        
                        context["density_mat"] = valss[0].tolist()
                        context["den_blocknone"] = ""
                        context["one_blocknone"] = ""                        
                    elif deriv == "radient":
                        rads = valss[0].tolist()
                        context["radient_mat"] = rads
                        context["rad_blocknone"] = ""
                        context["other_blocknone"] = ""                        
                    elif deriv == "laplacian":
                        laps = valss[0].tolist()
                        context["laplacian_mat"] = laps
                        context["lap_blocknone"] = ""
                        context["other_blocknone"] = ""                                            
                                                                                                                            
                
                print("Time taken to get get slice",datetime.datetime.now()-t1)
                t1 = datetime.datetime.now()
                ## Finally create the position dots if we want them
                pdots = settings_dic["posdots"]
                adots = settings_dic["atomdots"]                    
                print("Dots=",pdots,adots)
                context["zero_dotsX"] = []
                context["zero_dotsY"] = []
                context["posi_dotsX"] = []
                context["posi_dotsY"] = []                                        
                context["negi_dotsX"] = []
                context["negi_dotsY"] = []
                dots = []
                if pdots == "Y":
                    dots.append(central)
                    dots.append(linear)
                    dots.append(planar)
                if adots == "Y":
                    dots.append(cc)
                    dots.append(ll)
                    dots.append(pp)                    
                if len(dots) > 0:
                    ## Add points to a scatter plot
                    print(samples,width)
                    spc = space.SpaceTransform(central, linear, planar)
                    for dot in dots:
                        posD = spc.reverse_transformation(dot)                        
                        posDp = posD.get_point_pos(samples,width)                                                
                        #print("Dot",posDp.get_key())                                                
                        # The C value will be zero as it is on the plane - that is because these are the points we made the plane with
                        # The xy heatmap has been arranged so the x value is above so linear is upwards, so the y axis (ok a bit confusing.... should I change it)?
                        if posDp.A > 0 and posDp.A < samples and posDp.B > 0 and posDp.B < samples:
                            if abs(posDp.C) < 0.001:
                                context["zero_dotsX"].append(posDp.B)
                                context["zero_dotsY"].append(posDp.A)
                            elif posDp.C > 0:
                                context["posi_dotsX"].append(posDp.B)
                                context["posi_dotsY"].append(posDp.A)
                            elif posDp.C < 0:
                                context["negi_dotsX"].append(posDp.B)
                                context["negi_dotsY"].append(posDp.A)
                                                                                                                            
                print("Time taken to add dots",datetime.datetime.now()-t1)                
            except Exception as e:
                context["message"] = str(e)                    
                print(str(e))
                return 0            
        print("Time taken",datetime.datetime.now()-t1)
        print("Time taken",datetime.datetime.now()-ts)
        print("rendering...")        
        #################################################################################
        # it should be quicker now so loop and do it again

        return 0
    
    
    


###################################################################
## PROFILER ##
import cProfile
if include_timer:
    cProfile.run('runExploreScreen()')
else:
    # multiple times like refreshing when it already exists
    tms = 5
    ts = datetime.datetime.now()
    tss = []
    tse = []
    for tm in range(tms):
        tss.append(datetime.datetime.now())
        runExploreScreen()
        tse.append(datetime.datetime.now())

    print("each time")
    for t in range(tms):
        print(str(tse[t]-tss[t]))    
    print("Time taken to start check on file",datetime.datetime.now()-ts)


