"""
RSA 4/2/23


"""
from leuci_xyz import spacetransform as space
from leuci_xyz import crstransform as crs
from leuci_xyz import gridmaker as grid
from leuci_pol import interpolator as pol

#####################################################################
class MapFunctions(object):
    def __init__(self, pdb_code, mobj,pobj, interp_method):
        # PUBLIC INTERFACE
        self.pdb_code = pdb_code        
        self.mobj = mobj
        self.pobj = pobj
        self.interp = pol.create_interpolator(interp_method,self.mobj.values,self.mobj.F,self.mobj.M,self.mobj.S)
        self.interp_method = interp_method
        
    def get_slice(self,central, linear, planar, width, samples, interp_method, diff=0, log_level=0):
        # change interpolator if necessary
        if self.interp_method != interp_method:
            self.interp = pol.create_interpolator(interp_method,self.mobj.values,self.mobj.F,self.mobj.M,self.mobj.S)
        # data needed
        dim_order =  [int(self.mobj.map_header["01_NC"]),int(self.mobj.map_header["02_NR"]),int(self.mobj.map_header["03_NS"])]
        crs_starts =  [int(self.mobj.map_header["05_NCSTART"]),int(self.mobj.map_header["06_NRSTART"]),int(self.mobj.map_header["07_NSSTART"])]
        axis_sampling =  [int(self.mobj.map_header["08_NX"]),int(self.mobj.map_header["09_NY"]),int(self.mobj.map_header["10_NZ"])]
        map2crs =  [int(self.mobj.map_header["17_MAPC"])-1,int(self.mobj.map_header["18_MAPR"])-1,int(self.mobj.map_header["19_MAPS"])-1]
        cell_dims = [float(self.mobj.map_header["11_X_length"]),float(self.mobj.map_header["12_Y_length"]),float(self.mobj.map_header["13_Z_length"])]
        angles =  [float(self.mobj.map_header["14_Alpha"]),float(self.mobj.map_header["15_Beta"]),float(self.mobj.map_header["16_Gamma"])]
        #############        
        # objects needed
        spc = space.SpaceTransform(central, linear, planar)
        gm = grid.GridMaker()        
        #########                                        
        u_coords = gm.get_unit_grid(width,samples)        
        xyz_coords = spc.convert_coords(u_coords)
        crs_spc = crs.CrsTransform(dim_order, crs_starts, axis_sampling, map2crs, cell_dims, angles)        
        crs_coords = crs_spc.convert_coords_to_crs(xyz_coords)        
        vals = self.interp.get_val_slice(crs_coords,diff=diff)
        return vals
        



        
        
        
        