"""
RSA 10/3/23

This replicates the images submitted for my progress review, Birkbeck PhD March 2023, 2000 word status report
Using WSL the show() functions creates an html page on localhost

"""
################### USER INPUTS #######################
which_examples = [7,8] # could be 0,1,2,3,4 or 5
interp_method = "cubic"
width = 8           # in angstrom
samples = 100       # number of sample points along each axis to interpolate
########### A description of the examples #############
examples = []
#0, rings and light effect
examples.append(["0a","1yk4",2,-1,"(4.6,11.521,0.062)","(4.967,13.625,0.901)","(4.4801,13.4462,2.6594)",1.0,0.1])
#1ab, bond electrons
examples.append(["1a","1yk4",2,-1,"(4.021,16.6,0.386)","(3.335,15.784,1.472)","(3.493,16.727,-0.736)",1.0,0.3])
examples.append(["1b","1yk4",1,-1,"(4.021,16.6,0.386)","(3.335,15.784,1.472)","(3.493,16.727,-0.736)",0.8,0.8])
#2abcd, negative density in NOS switch
examples.append(["2a","3u7z",2,-1,"(37.849,30.718,1.938)","(40.218,30.095,1.15)","(24.16,38.119,-7.606)",1.0,0.8])
examples.append(["2b","3u7z",1,-1,"(37.849,30.718,1.938)","(40.218,30.095,1.15)","(24.16,38.119,-7.606)",0.8,0.8])

#4abcd, 1ejg distorted density for fo vs fc
examples.append(["4a","1ejg",2,-1,"(20.535,13.026,11.33)","(21.84,12.801,10.61)","(19.748,11.982,11.477)",1.0,0.8])
examples.append(["4b","1ejg",1,-1,"(20.535,13.026,11.33)","(21.84,12.801,10.61)","(19.748,11.982,11.477)",1.0,0.8])
examples.append(["4c","1ejg",1,0,"(20.535,13.026,11.33)","(21.84,12.801,10.61)","(19.748,11.982,11.477)",1.0,0.8])
examples.append(["4d","1ejg",0,1,"(20.535,13.026,11.33)","(21.84,12.801,10.61)","(19.748,11.982,11.477)",1.0,0.8])

##########################################################
## data
for which_example in which_examples:
  plotid = examples[which_example][0]
  pdb_code = examples[which_example][1]
  fo,fc = examples[which_example][2], examples[which_example][3]
  cc,cl,cp = examples[which_example][4],examples[which_example][5],examples[which_example][6]
  min_per, max_per = examples[which_example][7],examples[which_example][8]


  print(pdb_code,fo,fc,cc,cl,cp,min_per,max_per)
  ##########################################################
  ## Imports
  from pathlib import Path
  # phd/leuci-async/leuci-map
  CODEDIR = str(Path(__file__).resolve().parent.parent.parent)+ "/src/"
  import sys
  DATADIR = str(Path(__file__).resolve().parent.parent )+ "/data/"
  EG_DIR = str(Path(__file__).resolve().parent.parent )+ "/OfficialExamples/"
  sys.path.append(CODEDIR)

  from pathlib import Path
  # Ensure libraries are imported
  ###############################################################################
  from leuci_map import maploader as moad
  from leuci_map import mapfunctions as mfun
  from leuci_xyz import vectorthree as v3
  from leuci_xyz import spacetransform as space

  ##############################################################################

  # Create the map loader
  ###############################################################################
  ml = moad.MapLoader(pdb_code,directory=DATADIR)
  if not ml.exists():
      ml.download()
  ml.load()
  if ml.em_loaded:
    print("Loading values", pdb_code)
    ml.load_values()
    if ml.values_loaded:
      mf = mfun.MapFunctions(pdb_code,ml.mobj,ml.pobj,interp_method)
      print("Success loading", pdb_code)
    else:
      print("There is a problem loading",pdb_code)
  else:
    print("There is a problem loading",pdb_code)
  print("completed")
  ###############################################################################

  # Create the 3 coordinates for orientation
  ###############################################################################
  central = v3.VectorThree().from_coords(cc)
  linear = v3.VectorThree().from_coords(cl)
  planar = v3.VectorThree().from_coords(cp)
  print("completed")
  ###############################################################################

  # Create the 3 value slices
  ###############################################################################
  import datetime
  dt1 = datetime.datetime.now()
  print("Finding values and derivatives")
  print("density")
  vals = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=0,fo=fo,fc=fc)
  print("radient")
  rads = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=1,fo=fo,fc=fc)
  print("laplacian")
  laps = mf.get_slice(central,linear,planar,width,samples,interp_method,deriv=2,fo=fo,fc=fc)        
  dt2 = datetime.datetime.now()
  print("completed in", dt2-dt1)
  ###############################################################################
  ## Add points to a scatter plot
  spc = space.SpaceTransform(central, linear, planar)
  posC = spc.reverse_transformation(central)
  posL = spc.reverse_transformation(linear)
  posP = spc.reverse_transformation(planar)
  posCp = posC.get_point_pos(samples,width)
  posLp = posL.get_point_pos(samples,width)
  posPp = posP.get_point_pos(samples,width)
  scatterX = []
  scatterY = []
  if posCp.A > 0 and posCp.A < samples and posCp.B > 0 and posCp.B < samples:
    scatterX.append(posCp.B)
    scatterY.append(posCp.A)
  if posLp.A > 0 and posLp.A < samples and posLp.B > 0 and posLp.B < samples:
    scatterX.append(posLp.B)
    scatterY.append(posLp.A)
  if posPp.A > 0 and posPp.A < samples and posPp.B > 0 and posPp.B < samples:
    scatterX.append(posPp.B)
    scatterY.append(posPp.A)
  # The C value will be zero as it is on the plane - that is because these are the points we made the plane with

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

  data_scatter = go.Scatter(x=scatterX,y=scatterY,mode="markers",marker=dict(color="yellow",size=2),showlegend=False)
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
    data_lap = go.Heatmap(z=laps,colorscale=['Black','White'],showscale=False)
  else:
    data_den = go.Contour(z=vals,showscale=True, 
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
  fig.write_image(EG_DIR +"eg001_eg_" + plotid + ".jpg")
  print("#### Created image at", EG_DIR +"eg001_eg_" + plotid + ".jpg ####")