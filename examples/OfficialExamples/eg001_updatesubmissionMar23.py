"""
RSA 10/3/23

This replicates the images submitted for my progress review, Birkbeck PhD March 2023, 2000 word status report
Using WSL the show() functions creates an html page on localhost

"""
################### USER INPUTS #######################
which_examples = [11,15] # could be 0-14
interp_method = "bspline"
width = 8           # in angstrom
samples = 100       # number of sample points along each axis to interpolate
########### A description of the examples #############
examples = []
#0, rings and light effect
examples.append(["0a","1yk4",2,-1,"(4.6,11.521,0.062)","(4.967,13.625,0.901)","(4.4801,13.4462,2.6594)",1.0,0.1])     #0
#1ab, bond electrons
examples.append(["1a","1yk4",2,-1,"(4.021,16.6,0.386)","(3.335,15.784,1.472)","(3.493,16.727,-0.736)",1.0,0.3])       #1
examples.append(["1b","1yk4",1,-1,"(4.021,16.6,0.386)","(3.335,15.784,1.472)","(3.493,16.727,-0.736)",0.8,0.8])       #2
#2abcd, negative density in NOS switch
examples.append(["2a","3u7z",2,-1,"(37.849,30.718,1.938)","(40.218,30.095,1.15)","(24.16,38.119,-7.606)",1.0,0.8])    #3
examples.append(["2b","3u7z",1,-1,"(37.849,30.718,1.938)","(40.218,30.095,1.15)","(24.16,38.119,-7.606)",0.8,0.8])    #4
examples.append(["2c","3u7z",1,0,"(37.849,30.718,1.938)","(40.218,30.095,1.15)","(24.16,38.119,-7.606)",0.8,0.8])    #5
examples.append(["2d","3u7z",0,1,"(37.849,30.718,1.938)","(40.218,30.095,1.15)","(24.16,38.119,-7.606)",0.8,0.8])    #6
#3abcd, 1ejg distorted density for fo vs fc
examples.append(["3a","1ejg",2,-1,"(20.535,13.026,11.33)","(21.84,12.801,10.61)","(19.748,11.982,11.477)",1.0,0.8])   #7
examples.append(["3b","1ejg",1,-1,"(20.535,13.026,11.33)","(21.84,12.801,10.61)","(19.748,11.982,11.477)",1.0,0.8])   #8
examples.append(["3c","1ejg",1,0,"(20.535,13.026,11.33)","(21.84,12.801,10.61)","(19.748,11.982,11.477)",1.0,0.8])    #9
examples.append(["3d","1ejg",0,1,"(20.535,13.026,11.33)","(21.84,12.801,10.61)","(19.748,11.982,11.477)",1.0,0.8])    #10
#4abc occupancy
examples.append(["4a","6e6o",2,-1,"(24.334,6.71,5.327)","(23.503,6.977,4.467)","(24.0962,7.1878,6.7383)",1.0,0.8])   #11
examples.append(["4b","6e6o",1,-1,"(24.334,6.71,5.327)","(23.503,6.977,4.467)","(24.0962,7.1878,6.7383)",1.0,0.8])   #12
examples.append(["4c","6e6o",1,0,"(24.334,6.71,5.327)","(23.503,6.977,4.467)","(24.0962,7.1878,6.7383)",1.0,0.8])   #13
examples.append(["4d","6e6o",0,1,"(24.334,6.71,5.327)","(23.503,6.977,4.467)","(24.0962,7.1878,6.7383)",1.0,0.8])   #14
#5ab hydrogen bonds
examples.append(["5a","1r6j",2,-1,"(7.0004,-3.772,13.7901)","(4.6774,-4.601,15.2421)","(6.9464,-2.965,12.5851)",1.0,0.8])   #15
examples.append(["5b","1r6j",1,-1,"(7.0004,-3.772,13.7901)","(4.6774,-4.601,15.2421)","(6.9464,-2.965,12.5851)",1.0,0.8])   #16



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
  fo,fc = examples[which_example][2], examples[which_example][3]
  cc,cl,cp = examples[which_example][4],examples[which_example][5],examples[which_example][6]
  min_per, max_per = examples[which_example][7],examples[which_example][8]
  print(pdb_code,fo,fc,cc,cl,cp,min_per,max_per)
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
  ###############################################################################
  # Create the 3 value slices  
  print("Finding values and derivatives")
  print("density")
  vals = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=0,fo=fo,fc=fc)
  print("radient")
  rads = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=1,fo=fo,fc=fc)
  print("laplacian")
  laps = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=2,fo=fo,fc=fc)          
  ###############################################################################
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
  # Showing the plots in plotly
  # reference: https://plotly.com/python/contour-plots/
  ###############################################################################
  import plotly.graph_objs as go
  from plotly.subplots import make_subplots
  # mins and maxes for colors
  mind,maxd = 1000,-1000
  minl,maxl = 1000,-1000
  for i in range(len(laps)):
    for j in range(len(laps[0])):
      minl = min(laps[i][j],minl)
      maxl = max(laps[i][j],maxl)
      mind = min(vals[i][j],mind)
      maxd = max(vals[i][j],maxd)
  if maxd == mind:
    d0 = 0.5
  else:
    d0 = (0 - mind) / (maxd - mind)
  if minl == maxl:
    l0 = 0.5
  else:
    l0 = (0 - minl) / (maxl - minl)

  absmin = mind*min_per
  absmax = maxd*max_per
  for i in range(len(vals)):
    for j in range(len(vals[0])):        
      vals[i][j] = max(vals[i][j],absmin)
      vals[i][j] = min(vals[i][j],absmax)

  fig = make_subplots(rows=1, cols=3, subplot_titles=("Density", "Radient", "Laplacian"))

  data_scatter = go.Scatter(x=scatterX,y=scatterY,mode="markers",marker=dict(color="yellow",size=5),showlegend=False)
  if which_example == 0:
    data_den = go.Heatmap(z=vals,colorscale=['White','Black'],showscale=False)
    data_rad = go.Heatmap(z=rads,colorscale=['Black','White'],showscale=False)
    data_lap = go.Heatmap(z=laps,colorscale=['Black','White'],showscale=False)
  elif which_example == 2:
    data_den = go.Contour(z=vals,showscale=True, 
                          colorscale=[(0, "blue"), (d0, "grey"), (1, "crimson")],
                          line=dict(width=0.5,color="gray"))
    data_rad = go.Heatmap(z=rads,colorscale=['Black','White'],showscale=False)
    data_lap = go.Heatmap(z=laps,colorscale=['Black','White'],showscale=False)
  elif fo==1 and fc == -1:
    data_den = go.Contour(z=vals,showscale=True,
                          colorscale=[(0, "blue"), (d0, "snow"), (1, "crimson")],
                          line=dict(width=0.5,color="gray"))
    data_rad = go.Heatmap(z=rads,colorscale=['Black','White'],showscale=False)
    data_lap = go.Contour(z=laps,showscale=False, colorscale=[(0, "blue"), (l0, "snow"), (1, "crimson")])
  else:
    data_den = go.Contour(z=vals,showscale=False, 
                          colorscale=[(0, "grey"), (d0, "snow"), (d0+0.2, "cornflowerblue"),(0.9, "crimson"),(1.0, "rgb(100,0,0)")],
                          contours=dict(start=absmin,end=absmax,size=(absmax-absmin)/20),
                          line=dict(width=0.5,color="gray"))
    data_rad = go.Heatmap(z=rads,colorscale=['Black','White'],showscale=False)
    data_lap = go.Contour(z=laps,showscale=False, colorscale=[(0, "blue"), (l0, "snow"), (1, "crimson")])
  fig.add_trace(data_den,row=1,col=1)
  fig.add_trace(data_scatter,row=1,col=1)
  fig.add_trace(data_rad,row=1,col=2)
  fig.add_trace(data_scatter,row=1,col=2)
  fig.add_trace(data_lap,row=1,col=3)
  fig.add_trace(data_scatter,row=1,col=3)

  fig.update_xaxes(showticklabels=False) # hide all the xticks
  fig.update_yaxes(showticklabels=False) # hide all the xticks
  fig.update_yaxes(scaleanchor="x",scaleratio=1)
  fig.update_xaxes(scaleanchor="y",scaleratio=1)
  fig.write_image(EG_DIR +"eg001_eg_" + plotid + ".jpg",width=2000,height=1500)
  print("#### Created image at", EG_DIR +"eg001_eg_" + plotid + ".jpg ####")
  dt2 = datetime.datetime.now()
  print("completed in", dt2-dt1)