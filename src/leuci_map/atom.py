"""


"""

class Atom(object):
    def __init__(self, line):
        Line = line
        AtomNo = line.Substring(6, 5).Trim()
        AtomType = line.Substring(12, 4).Trim()
        occupant = line.Substring(16, 1).Trim()
        AA = line.Substring(17, 3).Trim()
        Chain = line.Substring(21, 1).Trim()
        ResidueNo = line.Substring(22, 5).Trim()
        insertion = line.Substring(26, 1).Trim()
        strX = line.Substring(30, 8).Trim()
        strY = line.Substring(38, 8).Trim()
        strZ = line.Substring(46, 8).Trim()
        X = Convert.ToDouble(strX)
        Y = Convert.ToDouble(strY)
        Z = Convert.ToDouble(strZ)
        Occupancy = Convert.ToDouble(line.Substring(54, 6).Trim())
        BFactor = Convert.ToDouble(line.Substring(60, 6).Trim())
        Element = line.Substring(66).Trim()
        Summary = ((((((Chain + ":") + ResidueNo) + ":") + AA) + ":") + AtomType)
    def distance(self, atm2):
        v1 = VectorThree(X, Y, Z)
        v2 = VectorThree(atm2.X, atm2.Y, atm2.Z)
        return v1.distance(v2)
    def coords(self):
        return VectorThree(X, Y, Z)
