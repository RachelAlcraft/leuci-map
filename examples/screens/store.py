
"""
RSA this is copied from leuci-web and should maybe be in leuci-map
Thread safe implementation of the Singleton Pattern
https://medium.com/analytics-vidhya/how-to-create-a-thread-safe-singleton-class-in-python-822e1170a7f6
See also: https://python-patterns.guide/gang-of-four/singleton/
 """

import threading
import datetime

class Store:
    _instance = None
    _lock = threading.Lock()
    strge_loader = {}
    strge_loader_dates = {}
    strge_interper = {}
    strge_interper_dates = {}
    
    def __new__(cls):
        if not cls._instance:  # This is the only difference
            with cls._lock:
                if not cls._instance:
                    cls._instance = super().__new__(cls)
        return cls._instance
                        
    ## LOADER ##
    def exists_loader(cls,pdb_code):
        return pdb_code in cls.strge_loader    
    def get_loader(cls,pdb_code):
        dt = cls.strge_loader_dates[pdb_code]
        lo = cls.strge_loader[pdb_code]
        return lo,dt    
    def add_loader(cls,pdb_code, map_obj):        
        cls.strge_loader[pdb_code] = map_obj
        cls.strge_loader_dates[pdb_code] = datetime.datetime.now()

    ## Interper ##
    def exists_interper(cls,pdb_code):
        return pdb_code in cls.strge_interper    
    def get_interper(cls,pdb_code):
        dt = cls.strge_interper_dates[pdb_code]
        lo = cls.strge_interper[pdb_code]
        return lo,dt    
    def add_interper(cls,pdb_code, map_obj):        
        cls.strge_interper[pdb_code] = map_obj
        cls.strge_interper_dates[pdb_code] = datetime.datetime.now()
            
    def clear(cls):
        cls.strge_loader.clear()
        cls.strge_loader_dates.clear()
        cls.strge_interper.clear()
        cls.strge_interper_dates.clear()

    def print_interpers(cls):
        ret_text = ""            
        for pdb_code,map_fun in cls.strge_interper.items():        
            ret_text += "\n" + pdb_code + "\t" + map_fun.mobj.header_as_string[:30] + "\n"                    
            if len(map_fun.mobj.values) > 0:
                dt = cls.strge_interper_dates[pdb_code]
                ret_text += "vals=" + str(len(map_fun.mobj.values)) + "\t\t"
                ret_text += str(dt) + "\n"
        return ret_text

                    
    