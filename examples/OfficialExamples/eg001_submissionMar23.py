"""
RSA 10/3/23

This replicates the images submitted for my progress review, Birkbeck PhD March 2023, 2000 word status report
Using WSL the show() functions creates an html page on localhost

"""
################### USER INPUTS #######################
which_examples = [2] # could be 0-5
interp_method = "bspline"
width = 8           # in angstrom
samples = 80       # number of sample points along each axis to interpolate
########### A description of the examples #############
examples = []
append = "-" + interp_method
#0, rings and light effect
examples.append(["Fig01-rings"+append,"1yk4",
                  ["(4.6,11.521,0.062)","(4.967,13.625,0.901)","(4.4801,13.4462,2.6594)"],
                  [("density",2,-1,1.0,0.1,(1,1),"BW")],
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
  cc,cl,cp = examples[which_example][2][0],examples[which_example][2][1],examples[which_example][2][2]
  plots = examples[which_example][3]  
  plot_config = examples[which_example][4]
  names = examples[which_example][5]
  print(pdb_code,cc,cl,cp,plots,plot_config,names)
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
    mf = mfun.MapFunctions(pdb_code,ml.mobj,ml.pobj,interp_method)      
    loadeds[pdb_code] = mf
  ###############################################################################
  # Create the 3 coordinates for orientation  
  central = v3.VectorThree().from_coords(cc)
  linear = v3.VectorThree().from_coords(cl)
  planar = v3.VectorThree().from_coords(cp)
  ## Add points to a scatter plot
  spc = space.SpaceTransform(central, linear, planar)
  posC = spc.reverse_transformation(central)
  posL = spc.reverse_transformation(linear)
  posP = spc.reverse_transformation(planar)
  posCp = posC.get_point_pos(samples,width)
  posLp = posL.get_point_pos(samples,width)
  posPp = posP.get_point_pos(samples,width)
  print("Central scatter=",posCp.get_key())
  print("Linear scatter=",posLp.get_key())
  print("Planar scatter=",posPp.get_key())
  scatterX = []
  scatterY = []
  # The C value will be zero as it is on the plane - that is because these are the points we made the plane with
  # The xy heatmap has been arranged so the x value is above so linear is upwards, so the y axis (ok a bit confusing.... should I change it)?
  if posCp.A > 0 and posCp.A < samples and posCp.B > 0 and posCp.B < samples:
    scatterX.append(posCp.B)
    scatterY.append(posCp.A)
  if posLp.A > 0 and posLp.A < samples and posLp.B > 0 and posLp.B < samples:
    scatterX.append(posLp.B)
    scatterY.append(posLp.A)
  if posPp.A > 0 and posPp.A < samples and posPp.B > 0 and posPp.B < samples:
    scatterX.append(posPp.B)
    scatterY.append(posPp.A)  

  ###############################################################################
  cols=plot_config[1]      
  if cols == 2:
    fig = make_subplots(rows=plot_config[0], cols=plot_config[1],subplot_titles=(names),horizontal_spacing=0.05,vertical_spacing=0.05,column_widths=[0.5,0.5])
  else:
    fig = make_subplots(rows=plot_config[0], cols=plot_config[1],subplot_titles=(names),horizontal_spacing=0.05,vertical_spacing=0.05)
  
  
  for deriv,fo,fc,min_per,max_per, plot_pos,hue in plots:
    print("Plot details=",deriv,fo,fc,min_per,max_per)    
    vals = [[]]
    if deriv == "density":    
      vals = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=0,fo=fo,fc=fc)
    elif deriv == "radient":
      vals = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=1,fo=fo,fc=fc)
    elif deriv == "laplacian":
      vals = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=2,fo=fo,fc=fc)
    ###############################################################################    
    # Showing the plots in plotly
    # reference: https://plotly.com/python/contour-plots/
    ###############################################################################    
    # mins and maxes for colors
    mind,maxd = 1000,-1000    
    for i in range(len(vals)):
      for j in range(len(vals[0])):        
        mind = min(vals[i][j],mind)
        maxd = max(vals[i][j],maxd)
    if maxd == mind:
      d0 = 0.5
    else:
      d0 = (0 - mind) / (maxd - mind)
    
    absmin = mind*min_per
    absmax = maxd*max_per
    #for i in range(len(vals)):
    #  for j in range(len(vals[0])):        
    #    vals[i][j] = max(vals[i][j],absmin)
    #    vals[i][j] = min(vals[i][j],absmax)
    
    data_scatter = go.Scatter(x=scatterX,y=scatterY,mode="markers",marker=dict(color="yellow",size=5),showlegend=False)
    
    if hue == "BW":
      data_vals = go.Heatmap(z=vals,colorscale=['White','Black'],showscale=False,zmin=absmin,zmax=absmax)
    elif hue == "RB":
      data_vals = go.Contour(z=vals,showscale=False,
                            colorscale=[(0, "rgb(100,0,0)"),(0.1, "crimson"), (d0, "silver"), (0.9, "blue"),(1, "navy")],
                            contours=dict(start=absmin,end=absmax,size=(absmax-absmin)/20),
                            line=dict(width=0.5,color="darkgray"),
                            zmin=absmin,zmax=absmax)
    else:
      data_vals = go.Contour(z=vals,showscale=False, 
                          colorscale=[(0, "grey"), (d0, "snow"), (d0+0.2, "cornflowerblue"),(0.9, "crimson"),(1.0, "rgb(100,0,0)")],
                          contours=dict(start=absmin,end=absmax,size=(absmax-absmin)/20),
                          line=dict(width=0.5,color="gray"),
                          zmin=absmin,zmax=absmax)
        
    fig.add_trace(data_vals,row=plot_pos[0],col=plot_pos[1])    
    fig.add_trace(data_scatter,row=plot_pos[0],col=plot_pos[1])    
  
  fig.update_xaxes(showticklabels=False) # hide all the xticks
  fig.update_yaxes(showticklabels=False) # hide all the xticks
  fig.update_yaxes(scaleanchor="x",scaleratio=1)    
  fig.update_xaxes(scaleanchor="y",scaleratio=1)
  rows, cols =plot_config[0],plot_config[1]
  wdth = 2000
  hight = int(wdth * rows/cols)
  fig.write_image(EG_DIR +"eg001_eg_" + plotid + ".jpg",width=wdth,height=hight)
  print("#### Created image at", EG_DIR +"eg001_eg_" + plotid + ".jpg ####")
  dt2 = datetime.datetime.now()
  print("completed in", dt2-dt1)