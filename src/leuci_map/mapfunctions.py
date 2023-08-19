"""
RSA 4/2/23


"""
from leuci_xyz import spacetransform as space
from leuci_xyz import crstransform as crs
from leuci_xyz import gridmaker as grid
from leuci_pol import interpolator as pol
from operator import itemgetter
import datetime
import math

#####################################################################
class MapFunctions(object):
    def __init__(self, pdb_code, mobj,pobj, interp_method,as_sd=0,fo=2,fc=-1,log_level=0):
        # PUBLIC INTERFACE
        self.pdb_code = pdb_code        
        self.mobj = mobj
        self.pobj = pobj
        self.interper_vals = [] #default top using only the main density, which for xray is 2Fo-1Fc
        self.as_sd = as_sd
        self.fo=fo
        self.fc=fc
        self.interp_method = interp_method.lower()
        
        if interp_method.lower()=="none":
            self.interper = None
        else:            
            self.make_interper_if_needed(interp_method,log_level,self.fo,self.fc)            
            #self.interper = pol.create_interpolator(interp_method,self.main_interper_vals,self.mobj.F,self.mobj.M,self.mobj.S,log_level=1,as_sd=self.as_sd)                
        self.crs_spc = None
        # data needed
        dim_order =  [int(self.mobj.map_header["01_NC"]),int(self.mobj.map_header["02_NR"]),int(self.mobj.map_header["03_NS"])]
        crs_starts =  [int(self.mobj.map_header["05_NCSTART"]),int(self.mobj.map_header["06_NRSTART"]),int(self.mobj.map_header["07_NSSTART"])]
        axis_sampling =  [int(self.mobj.map_header["08_NX"]),int(self.mobj.map_header["09_NY"]),int(self.mobj.map_header["10_NZ"])]
        map2crs =  [int(self.mobj.map_header["17_MAPC"])-1,int(self.mobj.map_header["18_MAPR"])-1,int(self.mobj.map_header["19_MAPS"])-1]
        cell_dims = [float(self.mobj.map_header["11_X_length"]),float(self.mobj.map_header["12_Y_length"]),float(self.mobj.map_header["13_Z_length"])]
        angles =  [float(self.mobj.map_header["14_Alpha"]),float(self.mobj.map_header["15_Beta"]),float(self.mobj.map_header["16_Gamma"])]
        self.crs_spc = crs.CrsTransform(dim_order, crs_starts, axis_sampling, map2crs, cell_dims, angles)        
             
    def make_interper_if_needed(self,interp_method,log_level,fo,fc):        
        vs,ds = 1,0
        calc_fofc = False
        create_interp = False
        if self.interper_vals == []:
            create_interp = True
            #if self.mobj.diff_values == [] or fo==2 and fc == -1:
            if self.mobj.diff_has == 0 or fo==2 and fc == -1:
                for i in range(len(self.mobj.values)):
                    self.interper_vals.append(self.mobj.values[i])
            else:
                calc_fofc = True
                
        elif self.mobj.diff_has != 0:
        #elif self.mobj.diff_values == []:
            if (self.interp_method != interp_method or self.fo != fo or self.fc != fc):
                calc_fofc = True
                create_interp = True
        if calc_fofc:
            self.interper_vals = []
            vs, ds = 0,0
            vs, ds = fo, -1 * fo                
            vs = vs + fc
            ds = ds + (-2 * fc)
            
            for i in range(len(self.mobj.values)):
                self.interper_vals.append(vs*self.mobj.values[i] + ds*self.mobj.diff_values[i])
        
        if create_interp:
            if log_level > 0:
                print("New interper, Fos=",fo,"Fcs=",fc,"mains=",vs,"diffs=",ds)
            self.interper = pol.create_interpolator(interp_method,self.interper_vals,(self.mobj.F,self.mobj.M,self.mobj.S),log_level=log_level,as_sd=self.as_sd)
        self.interp_method = interp_method
        self.fo = fo
        self.fc = fc
                    
    def max_min(self):
        return self.interper.min,self.interper.max

    def get_slices(self,central, linear, planar, width, samples, interp_method, derivs=[0], fo=2,fc=-1,log_level=0,degree=-1,ret_type="np"):
        vals = []
        for deriv in derivs:
            val = self.get_slice(central, linear, planar, width, samples, interp_method, deriv, fo,fc,log_level,degree,ret_type=ret_type)
            vals.append(val)
        return vals

    def get_slice(self,central, linear, planar, width, samples, interp_method, depth_samples=1, deriv=0, fo=2,fc=-1,log_level=0, ret_type="np"):
        # change interpolator if necessary        
        self.make_interper_if_needed(interp_method,log_level,fo,fc)        
        #############        
        # objects needed
        spc = space.SpaceTransform(central, linear, planar)        
        gm = grid.GridMaker()        
        #########                                        
        u_coords = gm.get_unit_grid(width,samples,depth_samples=depth_samples)
        xyz_coords = spc.convert_coords(u_coords)
        crs_coords = self.crs_spc.convert_coords_to_crs(xyz_coords)
        if depth_samples > 1:    
            vals = self.interper.get_val_slice3d(crs_coords,deriv=deriv,ret_type=ret_type)        
            return vals,xyz_coords
        else:
            vals = self.interper.get_val_slice(crs_coords,deriv=deriv,ret_type=ret_type)        
            return vals
    
    def get_slice_neighbours(self,central, linear, planar, width, samples,rnge,log_level=0):
        #############        
        # objects needed
        spc = space.SpaceTransform(central, linear, planar)
        gm = grid.GridMaker()        
        #########                                        
        u_coords = gm.get_unit_grid(width,samples)                
        xyz_coords = spc.convert_coords(u_coords)                
        in_scope = self.pobj.get_inscope_atoms(central,width+rnge[1],log_level)                
        a,b,c = xyz_coords.shape()
        naybs = []
        for i in range(a):
            row = []
            for j in range(b):
                coord = xyz_coords.get(i,j)                
                nayb = self.pobj.get_neighbours(coord,rnge,in_scope)
                row.append(nayb)
            naybs.append(row)
        return naybs
            
    def get_xyz(self,crs_coords):
        xyz_coords = self.crs_spc.crs_to_xyz(crs_coords)        
        return xyz_coords
    
    def get_crs(self,xyz_coords):
        crs_coords = self.crs_spc.xyz_to_crs(xyz_coords)        
        return crs_coords

    def get_map_projection(self, sliced, xmin=-1,xmax=-1,ymin=-1,ymax=-1):
        vals = self.interper.get_projection(sliced.lower(), xmin,xmax,ymin,ymax)        
        return vals
        
    def get_atoms_projection(self, interp_method,log_level=0):
        
        minx,miny,minz = 1000,1000,1000
        maxx,maxy,maxz = -1000,-1000,-1000

        atoms = self.pobj.get_atom_coords()
        self.make_interper_if_needed(interp_method,log_level,2,-1)
        crss = []        
        for atm in atoms:
            crs_coord = self.get_crs(atm)
            val = self.interper.get_value(crs_coord.A,crs_coord.B,crs_coord.C)
            crss.append((crs_coord,val))
        # now sort on the values
        sorted_crs = sorted(crss,key=itemgetter(1))
        xs,ys,zs,vs = [],[],[],[]
        for co,va in sorted_crs:
            xs.append(co.A)
            ys.append(co.B)
            zs.append(co.C)
            vs.append(va)
                        
            minx = math.floor(min(minx,float(co.A)))
            maxx = math.ceil(max(maxx,float(co.A)))
            miny = math.floor(min(miny,float(co.B)))
            maxy = math.ceil(max(maxy,float(co.B)))
            minz = math.floor(min(minz,float(co.C)))
            maxz = math.ceil(max(maxz,float(co.C)))                     
        
        return xs,ys,zs,vs,(minx,maxx),(miny,maxy),(minz,maxz)
    
    def get_map_cross_section(self, sliced,layer):
        vals = self.interper.get_cross_section(sliced.lower(),layer)        
        return vals

        



        
        
        
        