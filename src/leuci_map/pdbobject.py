"""
RSA 4/2/23
https://pynative.com/make-python-class-json-serializable/#:~:text=Use%20toJSON()%20Method%20to%20make%20class%20JSON%20serializable&text=So%20we%20don't%20need,Python%20Object%20to%20JSON%20string.

"""

import os
from os.path import exists
import json

class PdbObject(object):
    def __init__(self, pdb_code):
        # PUBLIC INTERFACE        
        self.pdb_code = pdb_code       
        self.lines = [] # consists of tuple of 
        #ATOM      1  N   GLY A 707   -0003.101   9.984  20.822  1.00  1.81    N  
        #         aid  atm  aa  ch  rid      x         y     z      occ   bf ele     
    
    def get_first_three(self):
        atm1, atm2, atm3 = {},{},{}
        count = 0
        for atm in self.lines:
            if atm["atm"] == "CA" and (atm["version"] == "A" or atm["version"] == ""):
                atm1 = atm
                count += 1
            elif atm["atm"] == "C" and (atm["version"] == "A" or atm["version"] == ""):
                atm2 = atm
                count += 1
            elif atm["atm"] == "O" and (atm["version"] == "A" or atm["version"] == ""):
                atm3 = atm
                count += 1
            if count == 3:
                return atm1,atm2,atm3
        return atm1,atm2,atm3

    def get_coords_key(self,key):
        atm = self.get_atm_key(key)
        if atm == {}:
            return "(0,0,0)"
        return f"({atm['x']},{atm['y']},{atm['z']})"

    def get_coords(self,atm):
        if atm == {}:
            return "(0,0,0)"
        return f"({atm['x']},{atm['y']},{atm['z']})"

    def get_atm_key(self,key):
        #A:709@C.A
        sx = key.split("@")
        s01 = sx[0].split(":")
        s23 = sx[1].split(".")
        ch,rid,am,ver = s01[0],s01[1],s23[0],s23[1]
        for atm in self.lines:        
            if atm["chain"] == ch and atm["rid"] == rid and atm["atm"] == am and atm["version"] == ver:        
                return atm
        return {}
    
    def get_key(self,atm):
        if atm == {}:
            return ""
        #A:709@C.A
        return f"{atm['chain']}:{atm['rid']}@{atm['atm']}.{atm['version']}"
    
    def get_next_key(self,key, offset=1):
        try:
            if len(key) < 4:
                return ""
            atm = self.get_atm_key(key)
            ridn = int(atm["rid"]) + offset
            return f"{atm['chain']}:{ridn}@{atm['atm']}.{atm['version']}"
        except:
            return ""
                        
    # if we add lines from a cif file I will do something different    
    def add_line_string(self,line):
        #ATOM    435  OE1AGLU A  23      12.418  13.330   0.541  0.58  7.57           O  
        if "ATOM" in line[:8] or "HETATM" in line[:8]:
            aid = line[8:13].strip()
            am = line[13:16].strip()
            ver = line[16:17].strip()
            aa = line[17:21].strip()
            ch = line[21:23].strip()
            rid = line[23:29].strip()
            x = line[29:38].strip()
            y = line[38:46].strip()
            z = line[46:54].strip()
            occ = line[54:60].strip()
            bf = line[60:66].strip()
            ele = line[66:].strip()
            atm = {}
            atm["rid"] = rid
            atm["aid"] = aid
            atm["aa"] = aa
            atm["atm"] = am
            atm["chain"] = ch
            if ver == "":
                ver = "A"
            atm["version"] = ver            
            atm["x"] = x
            atm["y"] = y
            atm["z"] = z
            atm["occupancy"] = occ
            atm["bfactor"] = bf
            atm["element"] = ele        
            self.lines.append(atm)

    def toJson(self):
        return json.dumps(self, default=lambda o: o.__dict__)
    
    def fromJson(self, jsndic):    
        self.pdb_code = jsndic["pdb_code"]
        self.lines = jsndic["lines"]
        
        
                                
    