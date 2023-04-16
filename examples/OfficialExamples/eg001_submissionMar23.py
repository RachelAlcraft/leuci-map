"""
RSA 10/3/23

This replicates the images submitted for my progress review, Birkbeck PhD March 2023, 2000 word status report
Using WSL the show() functions creates an html page on localhost

"""
################### USER INPUTS #######################
which_examples = [6] # could be 0-6
interp_method = "bspline"
width = 8           # in angstrom
samples = 100       # number of sample points along each axis to interpolate
########### A description of the examples #############
examples = []
append = "-" + interp_method
#0, rings and light effect
examples.append(["Fig01-rings"+append,"1yk4",
                  ["(4.6,11.521,0.062)","(4.967,13.625,0.901)","(4.4801,13.4462,2.6594)"],
                  [("density",2,-1,1.0,0.1,(1,1),"WB")],
                  (1,1),
                  ("density","")])     #0
#1ab, bond electrons
examples.append(["Fig02-bond-electrons"+append,"1yk4",
                  ["(4.021,16.6,0.386)","(3.335,15.784,1.472)","(3.493,16.727,-0.736)"],
                  [("density",2,-1,0.8,0.8,(1,1),"RGB"),("density",1,-1,0.4,0.4,(1,2),"RB"),
                  ("density",1,0,0.8,0.8,(2,1),"RGB"),("density",1,0,0.8,0.8,(2,2),"RGB")],
                  (2,2),
                  ("2FoFc","FoFc","Fo","Fc")])       #1
#2abcd, negative density in NOS switch
examples.append(["Fig03-NOS"+append,"3u7z",
                  ["(37.849,30.718,1.938)","(40.218,30.095,1.15)","(24.16,38.119,-7.606)"],
                  [("density",2,-1,0.8,0.5,(1,1),"RGB"),("density",1,-1,0.8,0.8,(1,2),"RB"),
                  ("density",1,0,0.8,0.5,(2,1),"RGB"),("density",0,1,0.8,0.5,(2,2),"RGB")],
                  (2,2),
                  ("2FoFc","FoFc","Fo","Fc")])       #2
#3abcd, 1ejg distorted density for fo vs fc
examples.append(["Fig04-distorted"+append,"1ejg",
                  ["(20.535,13.026,11.33)","(21.84,12.801,10.61)","(19.748,11.982,11.477)"],
                  [("density",2,-1,0.8,0.5,(1,1),"RGB"),("density",1,-1,0.8,0.8,(1,2),"RB"),
                  ("density",1,0,0.8,0.5,(2,1),"RGB"),("density",0,1,0.8,0.5,(2,2),"RGB")],
                  (2,2),
                  ("2FoFc","FoFc","Fo","Fc")])       #3
#4abc occupancy
examples.append(["Fig05-occupancy"+append,"6e6o",
                  ["(24.334,6.71,5.327)","(23.503,6.977,4.467)","(24.0962,7.1878,6.7383)"],
                  [("density",2,-1,1.0,0.8,(1,1),"RGB"),("density",1,-1,0.9,0.9,(1,2),"RB"),
                  ("density",1,0,1.0,0.8,(2,1),"RGB"),("density",0,1,1.0,0.8,(2,2),"RGB"),
                  ("radient",2,-1,1.0,1.0,(3,1),"BW"),("laplacian",2,-1,1.0,1.0,(3,2),"RB")],
                  (3,2),
                  ("2FoFc","FoFc","Fo","Fc","radient","laplacian")])       #4
#5ab hydrogen bonds
examples.append(["Fig06-hb"+append,"1r6j",
                  ["(7.0004,-3.772,13.7901)","(4.6774,-4.601,15.2421)","(6.9464,-2.965,12.5851)"],
                  [("density",2,-1,0.8,0.5,(1,1),"RGB"),("density",1,-1,0.8,0.8,(1,2),"RB"),
                  ("density",1,0,0.8,0.5,(2,1),"RGB"),("density",0,1,0.8,0.5,(2,2),"RGB")],
                  (2,2),
                  ("2FoFc","FoFc","Fo","Fc")])       #5

#7uly trial
examples.append(["Fig07-7uly"+append,"7uly",
                  ["(-6.172,-15.897,3.134)","(-6.447,-15.681,1.689)","(-5.583,-16.994,3.507)"],
                  [("density",2,-1,1.0,0.8,(1,1),"RGB"),("radient",2,-1,1.0,1.0,(1,2),"BW"),("laplacian",1,0,0.9,0.9,(1,3),"BR"),],
                  (1,3),
                  ("density","radient","laplacian")])       #6

########### Imports ################################
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


loadeds = {}
##########################################################
## data
for which_example in which_examples:
  dt1 = datetime.datetime.now()
  plotid = examples[which_example][0]
  pdb_code = examples[which_example][1]    
  cc,cl,cp = examples[which_example][2][0],examples[which_example][2][1],examples[which_example][2][2]
  plots = examples[which_example][3]  
  plot_config = examples[which_example][4]
  names = examples[which_example][5]
  print(pdb_code,cc,cl,cp,plots,plot_config,names)
  ##########################################################  
  ## Create the map loader  
  dt1 = datetime.datetime.now()
  if pdb_code in loadeds:
    print("Reusing loaded",pdb_code)
    mf = loadeds[pdb_code]
  else:
    print("Creating loaded",pdb_code)
    ml = moad.MapLoader(pdb_code,directory=DATADIR)
    if not ml.exists():
        ml.download()
    ml.load()
    if ml.em_loaded:
      print("Loading values", pdb_code)
      ml.load_values()
    if not ml.values_loaded:
      print("!!!! There is a problem loading",pdb_code)    
    mf = mfun.MapFunctions(pdb_code,ml.mobj,ml.pobj,interp_method)          
  ###############################################################################  
  # Create the 3 coordinates for orientation  
  central = v3.VectorThree().from_coords(cc)
  linear = v3.VectorThree().from_coords(cl)
  planar = v3.VectorThree().from_coords(cp)
  ################################################################################
  # Now create the MapPlotter
  filename = EG_DIR +"eg001_eg_" + plotid + ".jpg"
  mplot = mpl.MapPlotter(mf, filename,interp_method,samples,width,central,linear,planar,plot_config)
  ## Add each plot to the map    
  for deriv,fo,fc,min_per,max_per, plot_pos,hue in plots:
    print("Plot details=",deriv,fo,fc,min_per,max_per)
    mplot.add_plot_slice(deriv,fo,fc,min_per, max_per,hue,True,deriv,plot_pos)    
    
  mplot.make_plot_slices(log_level=1)
  print("#### Created image at", filename)
  dt2 = datetime.datetime.now()
  print("completed in", dt2-dt1)