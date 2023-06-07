"""
Created 21/5/2023 - Rachel Alcraft
Updates:

--------------------------------------
This script takes 
- a pdb file
- 3 atoms
- a width
- an interpolation method
And produces and html file with an inmage of the density using plotly

"""

#### USER INPUTS ###
pdb_code = "6eex"
central_atom = "A:707@C.A"
linear_atom = "A:707@CA.A"
planar_atom = "A:707@O.A"
interpolation = "bspline"
width = 8 #Angstrom
samples = 40
depth_samples = 10
file_output2d = "slice-2d.html"
file_output3d = "slice-3d.html"
#### ----------- ###




### CODE ###
from pathlib import Path
import sys
import matplotlib
CODEDIR = str(Path(__file__).resolve().parent.parent )+ "/src/"
sys.path.append(CODEDIR)
from leuci_xyz import vectorthree as v3
import leuci_map.mapobject as mobj
import leuci_map.maploader as moad
import leuci_map.mapfunctions as mfun
import leuci_map.mapsmanager as mman
#import leuci_map.mapplotter as mpl
import leuci_map.mapplothelp as mph

# find relative path for data
from pathlib import Path
DATADIR = str(Path(__file__).resolve().parent )+ "/data/"
RESDIR = str(Path(__file__).resolve().parent )+ "/results/" 
mman.MapsManager().set_dir(DATADIR)

# downloand/upload into memory the pdb and ccp4 data (0 means skip if not there, 1 means in this thread, 2 means in another thread and don't wait)
ml = mman.MapsManager().get_or_create(pdb_code,file=1,header=1,values=1)

mf = mfun.MapFunctions(pdb_code,ml.mobj,ml.pobj,interpolation)

cc = v3.VectorThree().from_coords(ml.pobj.get_coords_key(central_atom))
ll = v3.VectorThree().from_coords(ml.pobj.get_coords_key(linear_atom))
pp = v3.VectorThree().from_coords(ml.pobj.get_coords_key(planar_atom))

# 2d plot (s)
filename = RESDIR + file_output2d
vals = mf.get_slice(cc,ll,pp,width,samples,interpolation,deriv=0)
mplot = mph.MapPlotHelp(filename)
# add naybs
mfunc = mfun.MapFunctions(pdb_code,ml.mobj,ml.pobj, "linear") #the default method is linear
naybs = mfunc.get_slice_neighbours(cc,ll,pp,width,samples,[0,1],log_level=1)                        
mplot.make_plot_slice_2d(vals,min_percent=0.9,max_percent=0.9,title="ad my",samples=samples,width=width,points=[cc,ll,pp],naybs=naybs)
mplot = mph.MapPlotHelp(filename + "_no.html")
mplot.make_plot_slice_2d(vals,min_percent=0.9,max_percent=0.9,title="ad my",samples=samples,width=width,points=[cc,ll,pp])
# 3d plot
filename = RESDIR + file_output3d
mplot = mph.MapPlotHelp(filename)
vals,coords = mf.get_slice(cc,ll,pp,width,samples,interpolation,deriv=0,depth_samples=depth_samples)
mplot.make_plot_slice_3d(vals,min_percent=0.9,max_percent=0.9,title="3d my")

