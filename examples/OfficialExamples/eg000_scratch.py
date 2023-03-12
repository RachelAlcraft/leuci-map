"""
RSA 11/3/23

Scratch plots

"""
################### USER INPUTS #######################
if True:
  ## Imports
  from pathlib import Path
  # phd/leuci-async/leuci-map
  CODEDIR = str(Path(__file__).resolve().parent.parent.parent)+ "/src/"
  import sys
  DATADIR = str(Path(__file__).resolve().parent.parent )+ "/data/"
  EG_DIR = str(Path(__file__).resolve().parent.parent )+ "/OfficialExamples/"
  sys.path.append(CODEDIR)

  from pathlib import Path
  
  scatterX = []
  scatterY = []
  scatterX.append(0)
  scatterY.append(0)
  scatterX.append(5)
  scatterY.append(5)
  scatterX.append(1)
  scatterY.append(2)
  scatterX.append(1)
  scatterY.append(3)
  scatterX.append(1)
  scatterY.append(2.1)
  # The C value will be zero as it is on the plane - that is because these are the points we made the plane with

  # Showing the plots in plotly
  # reference: https://plotly.com/python/contour-plots/
  ###############################################################################
  import plotly.graph_objs as go
  from plotly.subplots import make_subplots
  # mins and maxes for colors
  
  fig = make_subplots(rows=1, cols=1, subplot_titles=("Scratch"))

  data_scatter = go.Scatter(x=scatterX,y=scatterY,mode="markers",marker=dict(color="red",size=5),showlegend=False)    
  fig.add_trace(data_scatter,row=1,col=1)  
  fig.update_xaxes(showticklabels=False) # hide all the xticks
  fig.update_yaxes(showticklabels=False) # hide all the xticks
  fig.update_yaxes(scaleanchor="x",scaleratio=1)
  fig.write_image(EG_DIR +"eg000_eg_scratch.jpg")
  print("#### Created image at", EG_DIR +"eg000_eg_scratch.jpg")