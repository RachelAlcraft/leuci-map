"""
RSA 4/2/23


"""

import os
from os.path import exists
import urllib.request
from Bio.PDB.MMCIFParser import MMCIFParser        
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBParser import PDBParser

import struct
import numpy as np

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

from . import mapobject as mobj
from . import pdbobject as pobj
from . import mapfunctions as mfun

class MapLoader(object):
    def __init__(self, pdb_code, directory="", cif=False):
        # PUBLIC INTERFACE
        self.mobj = mobj.MapObject(pdb_code)
        self.pobj = pobj.PdbObject(pdb_code)
        
        # Private data
        self._ccp4_binary = None
        self._diff_binary = None        
        self.pdb_loaded = False        
        self.em_loaded = False       
        self.values_loaded = False
        self.values_loading = False
        self.has_both = True
        self.mfunc = None
        self.mean = 0
        self.std = 1
        self.mean_diff = 0
        self.std_diff = 1
        # PRIVATE INTERFACE
        self._directory = directory        
        self._cif=cif
        if cif:
            self._filepath = f"{directory}{pdb_code}.cif"
            self.mobj.pdb_link = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_code}.cif"
        else:        
            self._filepath = f"{directory}{pdb_code}.pdb"
            self.mobj.pdb_link = f"https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{pdb_code}.ent"
        
        self._filepath_ccp4 = f"{self._directory}{self.mobj.pdb_code}.ccp4"
        self._filepath_diff = f"{self._directory}{self.mobj.pdb_code}_diff.ccp4"
                                        
    def exists(self):
        if self.exists_pdb():
            if self.exists_map():
                return True
        
        return False
    
    def exists_pdb(self):
        if exists(self._filepath) and exists(self._filepath+".done.txt"):
            return True
        else:
            return False
    
    def exists_map(self):
        self.load_pdb()
        if 'x-ray' in self.mobj.exp_method:
            if exists(self._filepath_ccp4) and exists(self._filepath_diff) and exists(self._filepath_ccp4+".done.txt"):
                #self.em_loaded = True
                return True        
            else:
                return False
        elif 'electron' in self.mobj.exp_method.lower():
            if exists(self._filepath_ccp4) and exists(self._filepath_ccp4+".done.txt"):
                #self.em_loaded = True
                return True
            else:
                return False
        else:
            #self.em_loaded = True
            return True # it doesn;t NOT exists anyway
    
    def success(self):
        if exists(self._filepath) and exists(self._filepath_ccp4):
            return True
        else:
            return False

    def download(self):                
        if not self.exists_pdb():
            if not self.already_started_pdb():
                self.write_starting_pdb()
                self.download_pdb()
        if not self.exists_map():
            if not self.already_started_map():
                self.write_starting_map()
                self.download_map()
        self.write_done()
    
    def already_started_pdb(self):
        return exists(self._filepath+".start.txt")    
    def already_started_map(self):
        return exists(self._filepath_ccp4+".start.txt")
    def write_starting_pdb(self):
        with open(self._filepath + ".start.txt","w") as fs_pdb:
            fs_pdb.write("start")            
    def write_starting_map(self):        
        with open(self._filepath_ccp4 + ".start.txt","w") as fs_map:
            fs_map.write("start")
    def write_done(self):
        if exists(self._filepath+".start.txt"):
            os.remove(self._filepath+".start.txt")
        if exists(self._filepath_ccp4+".start.txt"):
            os.remove(self._filepath_ccp4+".start.txt")
        if exists(self._filepath):
            with open(self._filepath + ".done.txt","w") as fw_pdb:
                fw_pdb.write("done")
        else:
            with open(self._filepath + ".failed.txt","w") as fw_pdb:
                fw_pdb.write("failed")
        if exists(self._filepath_ccp4):
            with open(self._filepath_ccp4 + ".done.txt","w") as fw_map:
                fw_map.write("done")
        else:
            with open(self._filepath_ccp4 + ".failed.txt","w") as fw_map:
                fw_map.write("failed")
                
    def download_pdb(self):
        self._fetch_pdbdata()

    def download_map(self):
        if not self.pdb_loaded:
            self.load_pdb()            
        if 'x-ray' in self.mobj.exp_method:
            self._fetch_maplink_xray()            
        elif 'electron' in self.mobj.exp_method:            
            print("Now downloading em matrix")
            self._fetch_maplink_em()
        #self.em_loaded = True
    
    def load(self):
        self.pdb_loaded = self.load_pdb()
        #if 'x-ray' in self.mobj.exp_method:            
        self.em_loaded = self.load_map()
            
            
        

    def load_pdb(self):
        if self._cif:
            structure = MMCIFParser().get_structure(self.mobj.pdb_code, self._filepath)
        else:
            structure = PDBParser(PERMISSIVE=True).get_structure(self.mobj.pdb_code, self._filepath)                
            found_em = False
            with open(self._filepath,"r") as fr:
                lines = fr.readlines()
                for line in lines:
                    self.pobj.add_line_string(line)
                    #REMARK 900 RELATED ID: EMD-6240   RELATED DB: EMDB                              
                    if "REMARK 900 RELATED ID:" in line and "EMD-" in line and not found_em:
                        found_em=True
                        info = line.split(" ")
                        self.mobj.em_code = info[4]
                        self.em_code = self.mobj.em_code
                        self.has_both = False
                        self.mobj.em_link = f"https://www.ebi.ac.uk/emdb/{self.em_code}"
                
        self._struc_dict = MMCIF2Dict(self._filepath)
        self.mobj.resolution = structure.header["resolution"]
        self.mobj.exp_method = structure.header["structure_method"]

        # with the structure create a pdbobject
        return True

    def load_map(self):        
        try:
            with open(self._filepath_ccp4, mode='rb') as file:
                self._ccp4_binary = file.read()
            if self.has_both:
                with open(self._filepath_diff, mode='rb') as file:
                    self._diff_binary = file.read()
            
            self._create_mapheader(self._ccp4_binary)            
            
            return True
        except Exception as e:
            print("Error loading map", str(e))
            return False
    
    def wait_for_load(self,log_level=0):
        import time
        count = 0
        while self.values_loading:
            time.sleep(0.1)
            if (log_level > -0):
                count += 1
                if count%30 == 0:
                    print("Waiting for load values", len(self.mobj.values), len(self.mobj.diff_values))
        return self.values_loaded

    def load_values(self, diff=True, adj_zero=True):
        try:            
            self.values_loading = True
            self._create_mapvalues(False,adj_zero=adj_zero)
            if diff and self.has_both:
                self._create_mapvalues(True,adj_zero=adj_zero)
            if not self.has_both:
                self.mobj.diff_values = []
            self.values_loading = False
            self.values_loaded = True
        except Exception as e:
            print(str(e))
            self.values_loaded = False
            self.values_loading = False


    def get_or_make_func(self, interp):
        if self.mfunc == None:
            self.mfunc = mfun.MapFunctions(self.mobj.pdb_code,self.mobj,self.pobj, interp) #the default method is linear
        return self.mfunc
    #################################################
    ############ PRIVATE INTERFACE ##################
    #################################################
    def _fetch_pdbdata(self):
        try:
            print(self.mobj.pdb_link, self._filepath)            
            urllib.request.urlretrieve(self.mobj.pdb_link, self._filepath)                                
        except:            
            return False
        return True
        
    def _fetch_maplink_xray(self):
        try:
            if not exists(self._filepath_ccp4):            
                urllib.request.urlretrieve(self.mobj.ccp4_link, self._filepath_ccp4)
            if not exists(self._filepath_diff):            
                urllib.request.urlretrieve(self.mobj.diff_link, self._filepath_diff)
            self.mobj.em_code = self.mobj.pdb_code
            self.em_code = self.mobj.pdb_code
        except Exception as e:
            print(e)

                        
    def _fetch_maplink_em(self):     
        """
        cif file
        EMDB EMD-6240 'associated EM volume' . 
        pdb file
        REMARK 900 RELATED ID: EMD-6240   RELATED DB: EMDB                              
        """
        self.mobj.em_code = self.em_code.split("-")[1]
        self.mobj.diff_link = "" # there is no difference density for cryo-em   
        self.mobj.ccp4_link = "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-" + self.mobj.em_code + "/map/emd_" + self.mobj.em_code + ".map.gz"        
        self.mobj.em_link = f"https://www.ebi.ac.uk/emdb/EMD-{self.em_code}"
        self.has_both = False
        # Download the files        
        filepath_gz = self._filepath_ccp4 + ".gz"
        if not exists(self._filepath_ccp4):
            # need to unzip file
            print("getting",self._filepath_ccp4)
            if not exists(filepath_gz):
                # need to download zipped file
                print("getting",self.mobj.ccp4_link,"to",filepath_gz)
                try:
                    urllib.request.urlretrieve(self.mobj.ccp4_link, filepath_gz)
                except Exception as e:
                    print(str(e))
            # Unzipping
            import gzip
            import shutil
            with gzip.open(filepath_gz, 'rb') as f_in:
                with open(self._filepath_ccp4, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                
        

    def _create_mapheader(self, ccp4_binary):
        num_labels = 0
        num_sym = 0        
        headers = [] #https://www.ccp4.ac.uk/html/maplib.html#description
        xheaders = [] 
        self.mobj.header_as_string = ""
        headers.append(["01_NC","int",4])           # of Columns    (fastest changing in map)
        headers.append(["02_NR","int",4])           # of Rows
        headers.append(["03_NS","int",4])           # of Sections   (slowest changing in map)        
        headers.append(["04_MODE","int",4])         # Data type   0 = signed bytes (from-128 lowest to 127 highest) 1 = Integer*2 2 = Image stored as Reals 3 = Complex Integer*2 4 = Complex Reals 5 == 0
        headers.append(["05_NCSTART","int",4])      # Number of first COLUMN  in map
        headers.append(["06_NRSTART","int",4])      # Number of first ROW     in map
        headers.append(["07_NSSTART","int",4])      # Number of first SECTION in map
        headers.append(["08_NX","int",4])           # Number of intervals along X
        headers.append(["09_NY","int",4])           # Number of intervals along Y
        headers.append(["10_NZ","int",4])           # Number of intervals along Z
        headers.append(["11_X_length","double",4])  # Cell Dimensions (Angstroms)
        headers.append(["12_Y_length","double",4])  #             "
        headers.append(["13_Z_length","double",4])  #             "
        headers.append(["14_Alpha","double",4])     # Cell Angles     (Degrees)
        headers.append(["15_Beta","double",4])      #             "
        headers.append(["16_Gamma","double",4])     #             "
        headers.append(["17_MAPC","int",4])         # Which axis corresponds to Cols.  (1,2,3 for X,Y,Z)
        headers.append(["18_MAPR","int",4])         # Which axis corresponds to Rows   (1,2,3 for X,Y,Z)
        headers.append(["19_MAPS","int",4])         # Which axis corresponds to Sects. (1,2,3 for X,Y,Z)
        headers.append(["20_AMIN","double",4])      # Minimum density value
        headers.append(["21_AMAX","double",4])      # Maximum density value
        headers.append(["22_AMEAN","double",4])     # Mean    density value    (Average)
        headers.append(["23_ISPG","int",4])         # Space group number
        headers.append(["24_NSYMBT","int",4])       # Number of bytes used for storing symmetry operators
        headers.append(["25_LSKFLG","int",4])       # Flag for skew transformation, =0 none, =1 if foll
        for i in range(26,35):
            headers.append([str(i) + "_SKWMAT","double",4])       # Flag for skew transformation, =0 none, =1 if foll
        for i in range(35,38):
            headers.append([str(i) + "_SKWTRN","double",4])       # Flag for skew transformation, =0 none, =1 if foll
        for i in range(38,53):
            headers.append(["X","int",4])       # Flag for skew transformation, =0 none, =1 if foll
        headers.append(["53_MAP","string",4])       # Character string 'MAP ' to identify file type
        headers.append(["54_MACHST","int",4])       # Machine stamp indicating the machine type
        headers.append(["55_ARMS","double",4])       # Rms deviation of map from mean density
        headers.append(["56_NLABL","int",4])       # Number of labels being used
                
        i=0
        for header, typ,inc  in headers:
            val = ""
            if not header == "X":
                if typ == "int":
                    val = int.from_bytes(ccp4_binary[i:i+inc], byteorder='little', signed=True)
                    self.mobj.map_header[header] = val
                elif typ == "double":
                    val = struct.unpack('f', ccp4_binary[i:i+inc])[0]
                    self.mobj.map_header[header] = val
                elif typ == "string":
                    val = ccp4_binary[i:i+inc].decode("utf-8") 
                    self.mobj.map_header[header] = val
                    
                if len(header) > 7:
                    self.mobj.header_as_string += header + "\t" + str(val) + "\n"
                else:
                    self.mobj.header_as_string += header + "\t\t" + str(val) + "\n"
                
                if header == "24_NSYMBT":
                    num_sym = int(val/80)                              
                if header == "56_NLABL":
                    num_labels = int(val)
                                                                
            i+=inc
                
        for s in range(0,num_labels):
            xheaders.append([str(s+1) + "_LABEL","string",80])       # 10  80 character text labels (ie. A4 format)
        for s in range(num_labels,10):
            xheaders.append(["X","string",80])       # 10  80 character text labels (ie. A4 format)            
        for s in range(0,num_sym):
            xheaders.append([str(s+1) + "_SYM","string",80])       # 10  80 character text labels (ie. A4 format)                        
            
        
        for header, typ,inc  in xheaders:
            if not header == "X":
                if typ == "int":
                    val = int.from_bytes(ccp4_binary[i:i+inc], byteorder='little', signed=True)
                    self.mobj.map_header[header] = val
                elif typ == "double":
                    val = struct.unpack('f', ccp4_binary[i:i+inc])[0]
                    self.mobj.map_header[header] = val
                elif typ == "string":
                    val = ccp4_binary[i:i+inc].decode("utf-8") 
                    self.mobj.map_header[header] = val
                    
                if len(header) > 7:
                    self.mobj.header_as_string += header + "\t" + str(val) + "\n"
                else:
                    self.mobj.header_as_string += header + "\t\t" + str(val) + "\n"
            i+=inc
        
                        
    def _create_mapvalues(self, diff, adj_zero=True):
        use_binary = self._ccp4_binary
        vals = []
        #vals = {} #zeros alternative
        if diff:
            use_binary = self._diff_binary
        
        Blength = self.mobj.map_header["01_NC"] * self.mobj.map_header["02_NR"] * self.mobj.map_header["03_NS"]
        Bstart = len(self._ccp4_binary) - (4 * Blength)

        self.mobj.F = self.mobj.map_header["01_NC"]
        self.mobj.M = self.mobj.map_header["02_NR"]
        self.mobj.S = self.mobj.map_header["03_NS"]

        if diff:
            self.mobj.diff_values = np.zeros((self.mobj.F,self.mobj.M,self.mobj.S))
            self.mobj.diff_has = self.mobj.F*self.mobj.M*self.mobj.S
        else:
            self.mobj.values = np.zeros((self.mobj.F,self.mobj.M,self.mobj.S))

        count = 0
        for s in range(0,self.mobj.S): #slow
            for m in range(0,self.mobj.M):#medium
                for f in range(0,self.mobj.F):#fast
                    strt = Bstart+(count*4)
                    val = struct.unpack('f', use_binary[strt:strt+4])[0]
                    count += 1
                    if diff:
                        self.mobj.diff_values[f,m,s] = val
                    else:
                        self.mobj.values[f,m,s] = val
                            
        #if adj_zero:
        #    if diff:
        #        zro = 
        #        self.mobj.diff_values = np.zeros((self.mobj.F,self.mobj.M,self.mobj.S))
        #    else:
        #        self.mobj.values = np.zeros((self.mobj.F,self.mobj.M,self.mobj.S))

        """
        Blength = self.mobj.map_header["01_NC"] * self.mobj.map_header["02_NR"] * self.mobj.map_header["03_NS"]
        Bstart = len(self._ccp4_binary) - (4 * Blength)
        for i in range(0,Blength):
            strt = Bstart+(i*4)
            val = struct.unpack('f', use_binary[strt:strt+4])[0]
            vals.append(val)
            #if abs(val) > 0.001:                
            #    vals[i] = val
        if diff:
            self.mobj.diff_values = vals
        else:
            self.mobj.values = vals
        """



    def _create_mapdata_em(self):
        ccp4_link = ""
        em_link = ""


    def cleanup(self):
        pass
