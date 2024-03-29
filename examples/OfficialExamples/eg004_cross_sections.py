"""
RSA 10/3/23

This shows projections

"""
################### USER INPUTS #######################
which_examples = [0] # could be 0-6
interp_method="linear"
########### A description of the examples #############
examples = []
append = ""

examples.append(["Fig01-cross-6eex"+append,"6eex",                  
                  [
                    ("xy",0,(1,1)),("yz",0,(1,2)),("zx",0,(1,3)),
                    ("xy",-1,(2,1)),("yz",-1,(2,2)),("zx",-1,(2,3)),
                    ],
                  (2,3),
                  ("xy","yz","zx","xy","yz","zx")])       #0

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
from leuci_xyz import vectorthree as v3
from leuci_xyz import spacetransform as space


loadeds = {}
##########################################################
## data
for which_example in which_examples:
  dt1 = datetime.datetime.now()
  plotid = examples[which_example][0]
  pdb_code = examples[which_example][1]      
  plots = examples[which_example][2]  
  plot_config = examples[which_example][3]
  names = examples[which_example][4]
  print(pdb_code,plots,plot_config,names)
  ##########################################################  
  ## Create the map loader  
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
    mf = mfun.MapFunctions(pdb_code,ml.mobj,ml.pobj,"linear")      
    loadeds[pdb_code] = mf
  ###############################################################################
  ###############################################################################
  cols=plot_config[1]      
  if cols == 2:
    fig = make_subplots(rows=plot_config[0], cols=plot_config[1],subplot_titles=(names),horizontal_spacing=0.05,vertical_spacing=0.05,column_widths=[0.5,0.5])
  else:
    fig = make_subplots(rows=plot_config[0], cols=plot_config[1],subplot_titles=(names),horizontal_spacing=0.05,vertical_spacing=0.05)
        
  for map_slice,layer,plot_pos in plots:    
    ################################################                
    if layer == -1: #then we want max
      if slice == "xy":
        layer = mfun.mobj.S
      elif slice == "yz":
        layer = mfun.mobj.F
      elif slice == "zx":
        layer = mfun.mobj.M
    print("Plot details=",map_slice,layer)    
    vals = mf.get_map_cross_section(map_slice,layer)
    data_vals = go.Heatmap(z=vals,colorscale=[(0, "snow"), (0.2, "cornflowerblue"),(0.9, "crimson"),(1.0, "rgb(100,0,0)")],showscale=False)               
    fig.add_trace(data_vals,row=plot_pos[0],col=plot_pos[1])
    ################################################
    

    
  fig.update_xaxes(showticklabels=False) # hide all the xticks
  fig.update_yaxes(showticklabels=False) # hide all the xticks
  fig.update_yaxes(scaleanchor="x",scaleratio=1)    
  fig.update_xaxes(scaleanchor="y",scaleratio=1)
  rows, cols =plot_config[0],plot_config[1]
  wdth = 2000
  hight = int(wdth * rows/cols)
  fig.write_image(EG_DIR +"eg004__eg_" + plotid + ".jpg",width=wdth,height=hight)
  print("#### Created image at", EG_DIR +"eg002_eg_" + plotid + ".jpg ####")
  dt2 = datetime.datetime.now()
  print("completed in", dt2-dt1)