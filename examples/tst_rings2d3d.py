
import sys
import os
from pathlib import Path
CODEDIR = str(Path(__file__).resolve().parent.parent )+ "/src/"
import sys
DATADIR = str(Path(__file__).resolve().parent )+ "/data/"
CODEDIR = str(Path(__file__).resolve().parent.parent )+ "/src/"
sys.path.append(CODEDIR)

pdb_code = '1yk4'
interpolation = 'bspline'
central,linear,planar = 'A:54@FE.A','A:28@CA.A','A:28@O.A'

from leuci_map import mapfunctions as mfun
from leuci_map import mapsmanager as mman
from leuci_map import mapplothelp as mph
from leuci_xyz import vectorthree as v3


print("Showing pdb map details", pdb_code)

interpolation = 'linear'
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

width = 6 #Angstrom
samples = 10
depth_samples = 10

cc = v3.VectorThree().from_coords(ml.pobj.get_coords_key(central))
ll = v3.VectorThree().from_coords(ml.pobj.get_coords_key(linear))
pp = v3.VectorThree().from_coords(ml.pobj.get_coords_key(planar))

cc = v3.VectorThree(4.6,11.521,0.062)
ll = v3.VectorThree(4.967,13.625,0.901)
pp = v3.VectorThree(4.4801,13.4462,2.6594)
print(cc.A,cc.B,cc.C)
print(ll.A,ll.B,ll.C)
print(pp.A,pp.B,pp.C)

naybs = []
if True:
    naybs = mf.get_slice_neighbours(cc,ll,pp,width,samples,[0,0.5],log_level=1)
vals2d = mf.get_slice(cc,ll,pp,width,samples,interpolation,deriv=0,ret_type="2d")
mplot = mph.MapPlotHelp("rings2d.html")
mplot.make_plot_slice_2d(vals2d,min_percent=0.5,max_percent=0.1,
samples=samples,width=width,points=[cc,ll,pp],title=pdb_code,hue="GRB",
plot_type="contour",
naybs=naybs)

#vals3dx,coords = mf.get_slice(cc,ll,pp,width,samples,interpolation,deriv=0,depth_samples=depth_samples,ret_type="3d")
#mplot3 = mph.MapPlotHelp("rings3d.html")
#mplot3.make_plot_slice_3d(vals3dx,min_percent=1, max_percent=0.05,title=pdb_code)