
import sys
import os
from pathlib import Path
CODEDIR = str(Path(__file__).resolve().parent.parent )+ "/src/"
import sys
DATADIR = str(Path(__file__).resolve().parent )+ "/data/"
CODEDIR = str(Path(__file__).resolve().parent.parent )+ "/src/"
sys.path.append(CODEDIR)

from leuci_map import mapfunctions as mfun
from leuci_map import mapsmanager as mman
from leuci_map import mapplothelp as mph
from leuci_xyz import vectorthree as v3

pdb_code = '1ejg'
interpolation = 'linear'
min_r,max_r = 0.0,0.1
levels = 10
hue = 'MASK'
width = 6 #Angstrom
samples = 50
depth_samples = 10
plot = 'contour'

if pdb_code == "1yk4":
    cc = v3.VectorThree(4.6,11.521,0.062)
    ll = v3.VectorThree(4.967,13.625,0.901)
    pp = v3.VectorThree(4.4801,13.4462,2.6594)
elif pdb_code == "1ejg":
    cc = v3.VectorThree(9.833,7.934,8.843)
    ll = v3.VectorThree(8.012,8.869,8.821)
    pp = v3.VectorThree(8.069,11.595,9.568)

print("Showing pdb map details", pdb_code)


ml = mman.MapsManager().get_or_create(pdb_code,file=1,header=1,values=1)
mf = mfun.MapFunctions(pdb_code,ml.mobj,ml.pobj,interpolation)

print("Map information...")
print('\tPdbCode= ',ml.mobj.pdb_code)
print('\tEmCode= ',ml.mobj.em_code)
print('\tPdbLink= ',ml.mobj.pdb_link)
print('\tEbiLink= ',ml.mobj.ebi_link)
#print(po.map_code)
print('\tResolution= ',ml.mobj.resolution)
print('\tExp Method= ',ml.mobj.exp_method)

print("Verifying values loaded...")
print('\t',len(ml.mobj.values),"=",ml.mobj.map_header["01_NC"] * ml.mobj.map_header["02_NR"] * ml.mobj.map_header["03_NS"])
print('\t',ml.mobj.values[0][0][0])    
print('\t',ml.mobj.values[len(ml.mobj.values)-1][0][0])  






print(cc.A,cc.B,cc.C)
print(ll.A,ll.B,ll.C)
print(pp.A,pp.B,pp.C)

naybs = []


if False:
    naybs = mf.get_slice_neighbours(cc,ll,pp,width,samples,[0,0.5],log_level=1)
vals2d = mf.get_slice(cc,ll,pp,width,samples,interpolation,deriv=0,ret_type="2d")
mplot = mph.MapPlotHelp("rings2d.html")
mplot.make_plot_slice_2d(vals2d,min_percent=min_r,max_percent=max_r,
samples=samples,width=width,points=[cc,ll,pp],title=pdb_code,hue=hue,
plot_type=plot,
naybs=naybs,
levels=levels)

vals3dx,coords = mf.get_slice(cc,ll,pp,width,samples,interpolation,deriv=0,depth_samples=depth_samples,ret_type="3d")
mplot3 = mph.MapPlotHelp("rings3d.html")
mplot3.make_plot_slice_3d(vals3dx,min_percent=min_r, max_percent=max_r,
                        title=pdb_code,transparency="low",hue=hue,
                        levels=levels)