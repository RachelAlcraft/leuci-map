"""
RSA 4/2/23

This loads and examines a map file and it's corresponding cif file given the pdb code
It will automatically check what kind of electron ddensity is available - xray or cryo em
"""

## Ensure code is importaed in path
from pathlib import Path
CODEDIR = str(Path(__file__).resolve().parent.parent )+ "/src/"
import sys
sys.path.append(CODEDIR)
from pathlib import Path
import leuci_xyz.vectorthree as v3
import leuci_xyz.spacetransform as sptr

########## INPUTS #################
central = v3.VectorThree(1,2,3)
linear = v3.VectorThree(2,2,2)
planar = v3.VectorThree(3,2,3)

########## EXAMPLE #################
def tranform_and_back(central, linear, planar):
    print("Transforming", central, linear, planar)
    
    # Test VectorThree
    vv = v3.v3_add(v3.VectorThree(0,1,2),v3.VectorThree(1,1,1))
    print("VectorThree test: central=", vv.A, vv.B, vv.C)

    # Test SpaceTransform goes back and forth
    st = sptr.SpaceTransform(central, linear, planar, log=True)
    check = st.apply_transformation(v3.VectorThree(0,0,0))
    print("SpaceTransform test, central returns to", check.A, check.B, check.C)

    # Test reverse goes back
    check_back = st.reverse_transformation(check)
    print("SpaceTransform test, central reverses to", check_back.A, check_back.B, check_back.C)

    # Test any old thing - stress    
    check_start = v3.VectorThree(15,12.6,-0.77)
    check3 = st.apply_transformation(check_start)
    check4 = st.reverse_transformation(check3)
    print("SpaceTransform test, from", check_start.A, check_start.B, check_start.C)
    print("SpaceTransform test, to", check3.A, check3.B, check3.C)
    print("SpaceTransform test, back", check4.A, check4.B, check4.C)



    
        

tranform_and_back(central, linear, planar)
