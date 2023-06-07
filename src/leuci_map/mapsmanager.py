
"""
RSA 16/4/23
Thread safe implementation of the Singleton Pattern
https://medium.com/analytics-vidhya/how-to-create-a-thread-safe-singleton-class-in-python-822e1170a7f6
See also: https://python-patterns.guide/gang-of-four/singleton/

This is intended to be the only access point to load maps
-----------------------------------------------------

 """

import threading #https://stackoverflow.com/questions/2905965/creating-threads-in-python
import datetime

from . import maploader as moad

class MapsManager:
    _instance = None
    _lock = threading.Lock()
    strge_container = {}
    DATADIR = ""
    CACHE = -1
        
    def __new__(cls):
        if not cls._instance:  # This is the only difference
            with cls._lock:
                if not cls._instance:
                    cls._instance = super().__new__(cls)
        return cls._instance                            
    
    def set_dir(cls,dir):
        cls.DATADIR = dir    
    def set_cache(cls,cache):
        cls.CACHE = cache
    
    ##  Map Store ##    
    def get_or_create(cls,pdb_code,file=1,header=1,values=1): #0 is no, 1 is in thread, 2 is new thread
        po = None
        if cls.exists_map(pdb_code):
            po = cls.get_map(pdb_code)            
        else:
            po = moad.MapLoader(pdb_code, directory=cls.DATADIR, cif=False)
            cls.strge_container[pdb_code] = po
        if not po.exists():
            if file == 1:
                po.download()
            elif file == 2:
                thread = threading.Thread(target=po.download,args=[])
                thread.start()            
        if not po.em_loaded:            
            if header == 1:
                po.load()
            elif header == 2:
                thread = threading.Thread(target=po.load,args=[])
                thread.start()            
        if not po.values_loaded:
            if values == 1:
                po.load()
            elif values == 2:
                thread = threading.Thread(target=po.load_values,args=[])
                thread.start()            
            po.load_values(diff=True)            
        return po
                                
    ############################################
    def exists_map(cls,pdb_code):
        return pdb_code in cls.strge_container    
    def get_map(cls,pdb_code):        
        lo = cls.strge_container[pdb_code]
        return lo
    def add_map(cls,pdb_code, map_obj):        
        cls.strge_container[pdb_code] = map_obj                    
    def clear(cls):
        cls.strge_container.clear()        
    def print_maps(cls):
        ret_text = ""            
        for pdb_code,po in cls.strge_container.items():        
            ret_text += "\n" + pdb_code + "\t" + po.mobj.header_as_string[:5] + "\t"                    
            if len(po.mobj.values) > 0:                
                ret_text += "vals=" + str(len(po.mobj.values)) + "\t\t"                
        return ret_text

                    
    