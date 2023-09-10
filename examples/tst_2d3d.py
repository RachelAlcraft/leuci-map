#! /usr/bin/env python

import sys
import os
from pathlib import Path
CODEDIR = str(Path(__file__).resolve().parent.parent )+ "/src/"
import sys
DATADIR = str(Path(__file__).resolve().parent )+ "/data/"
CODEDIR = str(Path(__file__).resolve().parent.parent )+ "/src/"
sys.path.append(CODEDIR)

from leuci_geo import pdbloader as pl
from leuci_geo import pdbgeometry as pg
from leuci_map import mapsmanager as mman
from leuci_map import mapfunctions as mfun
from leuci_xyz import vectorthree as v3
from leuci_map import mapplothelp as mph
import numpy as np


def slices(args): 
    pdb = args[1]
    PDBDIR = f"{args[2]}/"
    print(pdb,PDBDIR)

    mman.MapsManager().set_dir(PDBDIR)

    print("Loading",pdb)
    pla = pl.PdbLoader(pdb,PDBDIR,cif=False)    
    po = pla.load_pdb()

    gm = pg.GeometryMaker([po])
    geos = ['SG:{SG@1}[dis|<3.5]','SG:{(N),(O)}[dis|<3.5]','SG:{(N),(O)}:{SG@1}']
    df = gm.calculateGeometry(geos,log=0)

    print("----- filter -----")    
    df = df.loc[(df['occ_SG:{(N),(O)}:{SG@1}'] == 1) ]

    file_outputcsv = f"sg_csv_{pdb}.csv"
    df.to_csv(file_outputcsv)
    print("Saved to", file_outputcsv)

    coords_list = []
    i_count = 0
    df_len = len(df.index)
    atom_key = 'info_SG:{(N),(O)}:{SG@1}'
    for i, row in df.iterrows():
        pdb_code = row['pdb_code']
        atoms = row[atom_key]
        print(i_count, "/",df_len,i,atoms)
        i_count += 1
        clp = atoms.split("(")    
        cen = clp[1][:-1]
        lin = clp[2][:-1]
        pla = clp[3][:-1]
        cens = cen.split("|")
        lins = lin.split("|")
        plas = pla.split("|")
        #central_atom = "A:707@C.A"
        cen_str = cens[0]+":"+cens[2]+"@"+cens[3]+".A"
        lin_str = lins[0]+":"+lins[2]+"@"+lins[3]+".A"
        pla_str = plas[0]+":"+plas[2]+"@"+plas[3]+".A"
        
        coords_list.append((pdb_code,cen_str,lin_str,pla_str,atoms))
        print("Added",pdb_code,cen_str,lin_str,pla_str,atoms)

        print(coords_list)

        interpolation = "linear"    
        samples,width,depth_samples = 50,10,10        
        count = 1    
        for pdb_code,cen,lin,pla,atoms in coords_list:
            print(count, "/", len(coords_list),pdb_code,cen,lin,pla)
            #try:                        
            ml = mman.MapsManager().get_or_create(pdb_code,file=1,header=1,values=1)
            if not ml.success():
                print(pdb_code, "has not loaded ccp4 succesfully")
            else:
                mf = mfun.MapFunctions(pdb_code,ml.mobj,ml.pobj,interpolation,as_sd=2)
                cc = v3.VectorThree().from_coords(ml.pobj.get_coords_key(cen))
                ll = v3.VectorThree().from_coords(ml.pobj.get_coords_key(lin))
                pp = v3.VectorThree().from_coords(ml.pobj.get_coords_key(pla))            
                # 2d plot (s)
                filename2 = f"2d_matrices_{count}_{pdb_code}"                
                vals2d = mf.get_slice(cc,ll,pp,width,samples,interpolation,deriv=0,ret_type="2d")
                print(type(vals2d))
                np.save(filename2,vals2d)                                
                filename2 = f"2d_html_{count}_{pdb_code}.html"
                mplot = mph.MapPlotHelp(filename2)                
                mplot.make_plot_slice_2d(vals2d,min_percent=0.9,max_percent=0.9,samples=samples,width=width,points=[cc,ll,pp],title=pdb_code+":"+atoms)

                # 3d plot
                filename3 = f"3d_matrices_{count}_{pdb_code}"                        
                vals3d,coords = mf.get_slice(cc,ll,pp,width,samples,interpolation,deriv=0,depth_samples=depth_samples,ret_type="3d")
                np.save(filename3,vals3d)
                filename3 = f"3d_html_{count}_{pdb_code}.html"                
                mplot = mph.MapPlotHelp(filename3)                

                print(vals3d)



                mplot.make_plot_slice_3d(vals3d,min_percent=0.9, max_percent=0.9,title=pdb_code+":"+atoms)

#####################################################################################################                
slices(["","5i9s","/home/rachel/phd/leuci-async/leuci-flow/pdbdata"])         





