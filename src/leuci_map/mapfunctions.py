"""
RSA 4/2/23


"""

import struct
from leuci_xyz import spacetransform as space
from leuci_xyz import gridmaker as grid
from leuci_pol import nearest as near


class MapFunctions(object):
    def __init__(self, pdb_code, mobj):
        # PUBLIC INTERFACE
        self.pdb_code = pdb_code        
        self.mobj = mobj        

    def get_slice(self,central, linear, planar, width, samples):
        # objects needed
        spc = space.SpaceTransform(central, linear, planar)
        gm = grid.GridMaker()
        itrp = near.Nearest(self.mobj.values,self.mobj.F,self.mobj.M,self.mobj.S) 
        #########                                        
        u_coords = gm.get_unit_grid(width,samples)        
        xyz_coords = spc.get_coords(u_coords)
        vals = itrp.get_val_slice(xyz_coords)
        return vals
        



        
        
        
        