"""
RSA 17/03/23
This helper function makes it easier to plot for demonstrations and examples, paritclarly on the colab page


"""
from leuci_xyz import spacetransform as space
import plotly.graph_objs as go
from plotly.subplots import make_subplots

#####################################################################
class MapPlotter(object):
    def __init__(self, mfun, filename,interp_method,samples,width,central,linear,planar,plot_shape,degree=-1):
        # PUBLIC INTERFACE
        self.mfun = mfun
        self.filename = filename
        self.interp_method = interp_method
        self.samples = samples
        self.width=width
        self.central = central
        self.linear = linear
        self.planar = planar
        self.plot_shape = plot_shape
        self.plot_details = []
        self.sub_names= []
        
    def add_plot_slice(self,deriv,fo,fc,min_percent, max_percent,hue,dots,sub_name, plot_pos):
        self.plot_details.append((deriv,fo,fc,min_percent,max_percent,hue,dots,plot_pos))
        self.sub_names.append(sub_name)

    def make_plot_slices(self, log_level=0):
        # First create the dots for the potitions as a scatter plot
        spc = space.SpaceTransform(self.central, self.linear, self.planar)
        posC = spc.reverse_transformation(self.central)
        posL = spc.reverse_transformation(self.linear)
        posP = spc.reverse_transformation(self.planar)
        posCp = posC.get_point_pos(self.samples,self.width)
        posLp = posL.get_point_pos(self.samples,self.width)
        posPp = posP.get_point_pos(self.samples,self.width)
        if log_level > 0:
            print("Central scatter=",posCp.get_key())
            print("Linear scatter=",posLp.get_key())
            print("Planar scatter=",posPp.get_key())
        scatterX = []
        scatterY = []
        # The C value will be zero as it is on the plane - that is because these are the points we made the plane with
        # The xy heatmap has been arranged so the x value is above so linear is upwards, so the y axis (ok a bit confusing.... should I change it)?
        if posCp.A > 0 and posCp.A < self.samples and posCp.B > 0 and posCp.B < self.samples:
            scatterX.append(posCp.B)
            scatterY.append(posCp.A)
        if posLp.A > 0 and posLp.A < self.samples and posLp.B > 0 and posLp.B < self.samples:
            scatterX.append(posLp.B)
            scatterY.append(posLp.A)
        if posPp.A > 0 and posPp.A < self.samples and posPp.B > 0 and posPp.B < self.samples:
            scatterX.append(posPp.B)
            scatterY.append(posPp.A)  
  ###############################################################################
        cols=self.plot_shape[1]      
        if cols == 2:
            fig = make_subplots(rows=self.plot_shape[0], cols=self.plot_shape[1],subplot_titles=(self.sub_names),horizontal_spacing=0.05,vertical_spacing=0.05,column_widths=[0.5,0.5])
        else:
            fig = make_subplots(rows=self.plot_shape[0], cols=self.plot_shape[1],subplot_titles=(self.sub_names),horizontal_spacing=0.05,vertical_spacing=0.05)
    
        for deriv,fo,fc,min_per,max_per,hue,dots,plot_pos in self.plot_details:
            print("Plot details=",deriv,fo,fc,min_per,max_per)    
            vals = [[]]
            if deriv == "density":    
                vals = self.mfun.get_slice(self.central,self.linear,self.planar,self.width,self.samples,self.interp_method,deriv=0,fo=fo,fc=fc)
            elif deriv == "radient":
                vals = self.mfun.get_slice(self.central,self.linear,self.planar,self.width,self.samples,self.interp_method,deriv=1,fo=fo,fc=fc)
            elif deriv == "laplacian":
                vals = self.mfun.get_slice(self.central,self.linear,self.planar,self.width,self.samples,self.interp_method,deriv=2,fo=fo,fc=fc)
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
            data_scatter = go.Scatter(x=scatterX,y=scatterY,mode="markers",marker=dict(color="yellow",size=5),showlegend=False)            
            if hue == "BW":
                data_vals = go.Heatmap(z=vals,colorscale=['Black','Snow'],showscale=False,zmin=absmin,zmax=absmax)
            elif hue == "WB":
                data_vals = go.Heatmap(z=vals,colorscale=['Snow','Black'],showscale=False,zmin=absmin,zmax=absmax)
            elif hue == "RB":
                data_vals = go.Contour(z=vals,showscale=False,
                                    colorscale=[(0, "rgb(100,0,0)"),(0.1, "crimson"), (d0, "silver"), (0.9, "blue"),(1, "navy")],
                                    contours=dict(start=absmin,end=absmax,size=(absmax-absmin)/20),
                                    line=dict(width=0.5,color="darkgray"),
                                    zmin=absmin,zmax=absmax)
            elif hue == "BR":
                data_vals = go.Contour(z=vals,showscale=False,
                                    colorscale=[(0, "navy"),(0.1, "blue"), (d0, "silver"), (0.9, "crimson"),(1, "rgb(100,0,0)")],
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
            if dots:
                fig.add_trace(data_scatter,row=plot_pos[0],col=plot_pos[1])    
        
        fig.update_xaxes(showticklabels=False) # hide all the xticks
        fig.update_yaxes(showticklabels=False) # hide all the xticks
        fig.update_yaxes(scaleanchor="x",scaleratio=1)    
        fig.update_xaxes(scaleanchor="y",scaleratio=1)
        rows, cols =self.plot_shape[0],self.plot_shape[1]
        wdth = 2000
        hight = int(wdth * rows/cols)
        if self.filename == "SHOW":
            fig.show()
        elif ".html" in self.filename:
            fig.write_html(self.filename)
        else:
            fig.write_image(self.filename,width=wdth,height=hight)
        


                