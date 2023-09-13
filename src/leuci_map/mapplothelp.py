"""
RSA 17/03/23
This helper function makes it easier to plot for demonstrations and examples, paritclarly on the colab page


"""
from leuci_xyz import spacetransform as space
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import matplotlib

#####################################################################
class MapPlotHelp(object):
    def __init__(self,filename):
        # PUBLIC INTERFACE        
        self.filename = filename                
        self.plots_2d = []
        self.plots_3d = []                
            
    def add_plot_slice_2d(self,vals,tag):
        self.plots_2d.append((vals,tag))
    
    def add_plot_slice_3d(self,vals,tag):
        self.plots_3d.append((vals,tag))
                    
    def make_plot_slice_3d(self,vals,min_percent=1, max_percent=1,
                            hue="GBR",title="Leucippus Plot 3d",levels=20,
                            transparency="medium"):
        """
        Takes a mat3d object and plots it in plotly

        Input
        ---------
        vals : mat3d
        min_percent : 1
        max_percent : 1
        hue : GBR/GRB/BW/WB/RB/BR/GI/IG/MASK
        title : "Leucippud Plot 3d
        levels : 20
        """

        #https://plotly.com/python/3d-isosurface-plots/
        #turn data into scatter for iso_surface
        xs = []
        ys = []
        zs = []
        values = []

        #fig = make_subplots(rows=1, cols=1,subplot_titles=[title],horizontal_spacing=0.05,vertical_spacing=0.05)
        
        minv = 1000
        maxv = -1000

        a,b,c = vals.shape()
        for i in range(a):
            for j in range(b):
                for k in range(a):
                    val = 0
                    if k < c:
                        val = vals.get(i,j,k=k)
                        minv = min(minv,val)
                        maxv = max(maxv,val)
                    xs.append(i)
                    ys.append(j)
                    zs.append(k)
                    values.append(val)
                                    
        
        #if hue == "RINGS":
        #    minv,maxv =  -0.5,1.5
        #    min_percent,max_percent = 1,1
        
        #if hue != "MASK":
        absmin,absmax,d0,d1,d2 = self.__get_levels__(min_percent,max_percent,minv,maxv,hue=="MASK")
        #else:
        #absmin,absmax,d0,d1,d2 = minv,maxv,1,1,1
        
        
        
        colorscale=self.__get_colors__(hue,d0,d1,d2,min_percent,max_percent,transparency=transparency)
                
        fig= go.Figure(data=go.Isosurface(
        x=xs,
        y=ys,
        z=zs,
        value=values,        
        colorscale=colorscale,
        showscale=True,
        showlegend=False,        
        opacity=0.6,
        surface_count=levels,
        caps=dict(x_show=False, y_show=False,z_show=False),
        isomin=absmin,
        isomax=absmax,
        ),)
        
        fig.update_xaxes(showticklabels=False,visible=False,scaleanchor="x",scaleratio=1) # hide all the xticks
        fig.update_yaxes(showticklabels=False,visible=False,scaleanchor="x",scaleratio=1) # hide all the yticks
                                
        fig.update_layout(title=dict(text=title),
                            scene = dict(
                            xaxis = dict(visible=False),
                            yaxis = dict(visible=False),
                            zaxis =dict(visible=False)
                        )
    )

                        
        #print(values)
        if self.filename == "SHOW":
            fig.show()
        elif self.filename == "FIG":
            return fig
        elif ".html" in self.filename:
            fig.write_html(self.filename)
        else:
            fig.write_image(self.filename,width=2000,height=2000)

    def make_plot_slice_2d(self,vals2d,plot_type="contour",points=[], 
                            naybs=[],min_percent=1, max_percent=1,
                            hue="GBR",levels=20,title="Leucippus Plot 2d",samples=-1,width=-1,
                            transparency="no"):
        """
        Takes a mat2d object and plots it in plotly

        Input
        ---------
        vals : mat2d
        plottype : "countour" or heatmap
        points : [] a selection of points to add ontop of plot
        naybs : [] neighbours matched per value points
        min_percent : 1
        max_percent : 1
        hue : GBR/GRB/BW/WB/RB/BR/GI/IG/MASK
        title : "Leucippud Plot 2d
        levels : 20
        samples :-1 needed to place the points in the right place
        width : -1 needed to place the points in the right place
        """        
        #https://plotly.com/python/3d-isosurface-plots/
        #turn data into scatter for iso_surface        
        vals = vals2d.tolist()
        #print(type(vals))
        
        fig = make_subplots(rows=1, cols=1,subplot_titles=[title],horizontal_spacing=0.05,vertical_spacing=0.05)

        mind,maxd = 1000,-1000    
        for i in range(len(vals)):
            for j in range(len(vals[0])):        
                mind = min(vals[i][j],mind)
                maxd = max(vals[i][j],maxd)        
        
        #if hue == "RINGS":
        #    mind,maxd = -0.5,1.5
        #    min_percent,max_percent = 1,1
        
        #if hue != "MASK":
        absmin,absmax,d0,d1,d2 = self.__get_levels__(min_percent,max_percent,mind,maxd,hue=="MASK")        
        #else:
        #    absmin,absmax,d0,d1,d2 = mind,maxd,1,1,1
        colorscale=self.__get_colors__(hue,d0,d1,d2,min_percent,max_percent,transparency=transparency,)

        if len(naybs) > 0:
            if plot_type == "contour":
                data_vals = go.Contour(z=vals,showscale=False, 
                                colorscale=colorscale,
                                contours=dict(start=absmin,end=absmax,size=(absmax-absmin)/levels),
                                text=naybs,   
                                hovertemplate='......%{z:.4f}<br>%{text}',
                                line=dict(width=0.5,color="gray"),
                                zmin=absmin,zmax=absmax,name='')
            elif plot_type == "heatmap":
                data_vals = go.Heatmap(z=vals,showscale=False, 
                                colorscale=colorscale,                                
                                text=naybs,   
                                hovertemplate='......%{z:.4f}<br>%{text}',                                
                                zmin=absmin,zmax=absmax,name='')
        else:
            if plot_type == "contour":
                data_vals = go.Contour(z=vals,showscale=True, 
                                colorscale=colorscale,
                                contours=dict(start=absmin,end=absmax,size=(absmax-absmin)/levels),
                                hovertemplate='......%{z:.4f}',
                                line=dict(width=0.5,color="gray"),
                                zmin=absmin,zmax=absmax,name='')
            elif plot_type == "heatmap":
                data_vals = go.Heatmap(z=vals,showscale=True, 
                                colorscale=colorscale,                                
                                hovertemplate='......%{z:.4f}',                                
                                zmin=absmin,zmax=absmax,name='')

                        
        fig.add_trace(data_vals,row=1,col=1)
        if len(points) == 3:
            data_scatter = self.add_points(points,samples,width)
            fig.add_trace(data_scatter,row=1,col=1)    
        fig.update_xaxes(showticklabels=False,visible=False) # hide all the xticks
        fig.update_yaxes(showticklabels=False,visible=False) # hide all the xticks
        fig.update_yaxes(scaleanchor="x",scaleratio=1)    
        fig.update_xaxes(scaleanchor="y",scaleratio=1)
                        
        #print(values)
        if self.filename == "SHOW":
            fig.show()
        elif self.filename == "FIG":
            return fig
        elif ".html" in self.filename:
            fig.write_html(self.filename)        
        else:
            fig.write_image(self.filename,width=2000,height=2000)

    def add_points(self, points,samples,width,log_level=0):
        # First create the dots for the potitions as a scatter plot
        spc = space.SpaceTransform(points[0], points[1], points[2])
        posC = spc.reverse_transformation(points[0])
        posL = spc.reverse_transformation(points[1])
        posP = spc.reverse_transformation(points[2])
        posCp = posC.get_point_pos(samples,width)
        posLp = posL.get_point_pos(samples,width)
        posPp = posP.get_point_pos(samples,width)
        if log_level > 0:
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
                    
        data_scatter = go.Scatter(x=scatterX,y=scatterY,mode="markers",marker=dict(color="yellow",size=5),showlegend=False,hoverinfo='skip',hovertemplate='',name='')            

        return data_scatter

    def __get_levels__(self,min_rr, max_rr, min_vv, max_vv,absolute):        
        
        if absolute:
            mn = min(min_rr,max_rr)
            mx = max(min_rr,max_rr)
            absmin = (mn*(max_vv-min_vv)) + min_vv
            absmax = (mx*(max_vv-min_vv)) + min_vv
        else:
            absmin = min_vv*min_rr        
            absmax = max_vv*max_rr          
            
        if absmin == absmax:
            d0 = 0.5
        elif absmin > 0:
            d0 = 0
        else:
            d0 = (0 - absmin) / (absmax - absmin)                            
        d1 = d0 + ((1-d0)/3)
        d2 = d0 + (2*(1-d0)/3)                        
        return absmin,absmax,d0,d1,d2
    
    def __get_colors__(self,hue,d0,d1,d2,min_r,max_r,transparency):
        t0,t1,t2,t3,t4,tx,tl,tm,tu = 1,1,1,1,1,1,1,1,1
        if transparency == "low":
            t0,t1,t2,t3,t4 = 0.5,0.5,0.5,0.5,0.5
            tx = 0.5
            tl,tm,tu = 0,0.5,1
        elif transparency == "medium":
            t0,t1,t2,t3,t4 = 0.4,0.1,0.4,0.5,0.4
            tx = 0.5
            tl,tm,tu = 0.1,0.5,0.4
        elif transparency == "high":
            t0,t1,t2,t3,t4 = 0,0.7,0.5,0.2,0
            tx = 0.5
            tl,tm,tu = 0.1,0.5,0.1                       
        if hue == 'GRB':
            t0,t1,t2,t3,t4 = t0,t1,t4,t3,t2
        elif hue == "RB" or hue == 'GI':
            tl,tm,tu = tu,tm,tl

        grey = f"rgba(119,136,153,{t0})"
        snow = f"rgba(240,248,255,{t1})"
        cornflower = f"rgba(100,149,237,{t2})"
        crimson = f"rgba(220,20,60,{t3})"
        darkred = f"rgba(100,0,0,{t4})"
        navy = f"rgba(0,0,128,{t4})"

        black = f"rgba(0,0,0,{tx})"
        ghost = f"rgba(248,248,255,{tx})"

        midnight = f"rgba(25,25,112,{tl})"
        maroon = f"rgba(128,0,0,{tu})"

        indigo = f"rgba(75,0,130,{tl})"        
        sea = f"rgba(32,178,170,{tu})"

        slate = f"rgba(40,79,79,{0.7})"        
        fire = f"rgba(178,34,34,{0.7})"  
        alice = f"rgba(240,248,255,{0.7})"  
        rose = f"rgba(255,228,225,{0.7})" 



        
        if hue == "GBR":            
            return [(0, grey), (d0, snow), (d1, cornflower),(d2, crimson),(1.0, darkred)]
        elif hue == "GRB":            
            return [(0, grey), (d0, snow), (d1, crimson),(d2, cornflower),(1.0, navy)]                
        elif hue == "BW":
            return [(0, black),(1.0, ghost)]
        elif hue == "WB":
            return [(0, ghost),(1.0, black)]
        elif hue == "RB":
            if d0 <= 0:
                return [(0,maroon),(0.5,ghost),(1.0,midnight)]
            else:
                return [(0,maroon),(d0,ghost),(1.0,midnight)]
        elif hue == "BR":
            if d0 <= 0:
                return [(0,midnight),(0.5,ghost),(1.0,maroon)]
            else:
                return [(0,midnight),(d0,ghost),(1.0,maroon)]
        elif hue == "IG":
            if d0 <= 0:
                return [(0,indigo),(0.5,ghost),(1.0,sea)]
            else:
                return [(0,indigo),(d0,ghost),(1.0,sea)]
        elif hue == "GI":
            if d0 <= 0:
                return [(0,sea),(0.5,ghost),(1.0,indigo)]
            else:
                return [(0,sea),(d0,ghost),(1.0,indigo)]
        elif hue == "MASK":
            m2 = min(min_r,max_r)
            m3 = max(min_r,max_r)
            m1 = max(0,min(min_r,max_r)-0.01)
            m4 = min(1,max(min_r,max_r)+0.01)
            #return [(0,shell),(m1,shell),(m2,slate),(m3,fire),(m4,shell),(1.0,shell)]
            return [(0,alice),(0.1,alice),(0.2,slate),(0.8,fire),(0.9,rose),(1.0,rose)]                        
            #return [(0,shell),(0.5,slate),(1.0,shell)]
            
            

      

        


        
    

            
    


                    