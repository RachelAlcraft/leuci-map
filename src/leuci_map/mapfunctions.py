"""
RSA 4/2/23


"""
from leuci_xyz import spacetransform as space
from leuci_xyz import crstransform as crs
from leuci_xyz import gridmaker as grid
from leuci_pol import interpolator as pol
from operator import itemgetter

#####################################################################
class MapFunctions(object):
    def __init__(self, pdb_code, mobj,pobj, interp_method,degree=-1):
        # PUBLIC INTERFACE
        self.pdb_code = pdb_code        
        self.mobj = mobj
        self.pobj = pobj
        self.main_interper_vals = [] #default top using only the main density, which for xray is 2Fo-1Fc
        if interp_method.lower()=="none":
            self.interper = None
        else:
            for i in range(len(self.mobj.values)):
                self.main_interper_vals.append(self.mobj.values[i])
            self.interper = pol.create_interpolator(interp_method,self.main_interper_vals,self.mobj.F,self.mobj.M,self.mobj.S,log_level=1,degree=degree)
        self.interp_method = interp_method
        self.degree = degree
        self.fo=2
        self.fc=-1
        self.crs_spc = None
        # data needed
        dim_order =  [int(self.mobj.map_header["01_NC"]),int(self.mobj.map_header["02_NR"]),int(self.mobj.map_header["03_NS"])]
        crs_starts =  [int(self.mobj.map_header["05_NCSTART"]),int(self.mobj.map_header["06_NRSTART"]),int(self.mobj.map_header["07_NSSTART"])]
        axis_sampling =  [int(self.mobj.map_header["08_NX"]),int(self.mobj.map_header["09_NY"]),int(self.mobj.map_header["10_NZ"])]
        map2crs =  [int(self.mobj.map_header["17_MAPC"])-1,int(self.mobj.map_header["18_MAPR"])-1,int(self.mobj.map_header["19_MAPS"])-1]
        cell_dims = [float(self.mobj.map_header["11_X_length"]),float(self.mobj.map_header["12_Y_length"]),float(self.mobj.map_header["13_Z_length"])]
        angles =  [float(self.mobj.map_header["14_Alpha"]),float(self.mobj.map_header["15_Beta"]),float(self.mobj.map_header["16_Gamma"])]
        self.crs_spc = crs.CrsTransform(dim_order, crs_starts, axis_sampling, map2crs, cell_dims, angles)        
        
    
    def make_interper_if_needed(self,interp_method,log_level,fo,fc,degree):
        change_degree = (degree != -1) and (degree != self.degree)
        if (self.interp_method != interp_method or change_degree) and self.mobj.diff_values == []:#then the method has changed but it is em with no diffs
            self.interper = pol.create_interpolator(interp_method,self.main_interper_vals,self.mobj.F,self.mobj.M,self.mobj.S,degree=degree,log_level=log_level)
        elif (self.interp_method != interp_method or change_degree) or (self.fo != fo or self.fc != fc and self.mobj.diff_values != []):            
            if self.mobj.diff_values != []: #we can't change the fo and fc if thee is no diff
                interper_vals = []
                vs, ds = 0,0
                vs, ds = fo, -1 * fo                
                vs = vs + fc
                ds = ds + (-2 * fc)
                if log_level > 0:
                    print("New interper, Fos=",fo,"Fcs=",fc,"mains=",vs,"diffs=",ds)
                for i in range(len(self.mobj.values)):
                    interper_vals.append(vs*self.mobj.values[i] + ds*self.mobj.diff_values[i])
                self.interper = pol.create_interpolator(interp_method,interper_vals,self.mobj.F,self.mobj.M,self.mobj.S,log_level=log_level)
        self.interp_method = interp_method
        self.fo = fo
        self.fc = fc
        if change_degree:
            self.degree = degree
    
    def get_slice(self,central, linear, planar, width, samples, interp_method, deriv=0, fo=2,fc=-1,log_level=0,degree=-1):
        # change interpolator if necessary
        self.make_interper_if_needed(interp_method,log_level,fo,fc,degree)                
        #############        
        # objects needed
        spc = space.SpaceTransform(central, linear, planar)
        gm = grid.GridMaker()        
        #########                                        
        u_coords = gm.get_unit_grid(width,samples)        
        xyz_coords = spc.convert_coords(u_coords)        
        crs_coords = self.crs_spc.convert_coords_to_crs(xyz_coords)        
        vals = self.interper.get_val_slice(crs_coords,deriv=deriv)
        return vals
    
    def get_xyz(self,crs_coords):
        xyz_coords = self.crs_spc.crs_to_xyz(crs_coords)        
        return xyz_coords
    
    def get_crs(self,xyz_coords):
        crs_coords = self.crs_spc.xyz_to_crs(xyz_coords)        
        return crs_coords

    def get_map_projection(self, sliced):
        vals = self.interper.get_projection(sliced.lower())        
        return vals
    
    def get_atoms_projection(self, interp_method,degree=-1,log_level=0):
        atoms = self.pobj.get_atom_coords()
        self.make_interper_if_needed(interp_method,log_level,2,-1,degree)
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
        return xs,ys,zs,vs
    
    def get_map_cross_section(self, sliced,layer):
        vals = self.interper.get_cross_section(sliced.lower(),layer)        
        return vals

        



        
        
        
        