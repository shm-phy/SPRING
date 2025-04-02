
import argparse
import numpy as np


def FindLatVec(a, b, c, angl_bc, angl_ac, angl_ab) :

    L1 = np.zeros(3, dtype=np.double)
    L2 = np.zeros(3, dtype=np.double)
    L3 = np.zeros(3, dtype=np.double)

    L1[0] = a

    L2[0] = b * np.cos( angl_ab )
    L2[1] = b * np.sin( angl_ab )
    
    L3[0] = c * np.cos( angl_ac )
    L3[1] = c * ( (np.cos(angl_bc) - np.cos(angl_ab) * np.cos(angl_ac))  / np.sin(angl_ab) )
    L3[2] = np.sqrt( c**2 - L3[0]**2 - L3[1]**2 )

    return L1, L2, L3


if __name__ == "__main__" :

    #Run: python3 LatticeVectors.py -a 3.952048801 -b 3.952048801 -c 3.952048801 -alpha 60 -beta 60 -gama 60

    parser = argparse.ArgumentParser()

    parser.add_argument("-a", type=float, default=1.0, help="Lattice constant a in Anstrom")
    parser.add_argument("-b", type=float, default=1.0, help="Lattice constant b in Anstrom")
    parser.add_argument("-c", type=float, default=1.0, help="Lattice constant c in Anstrom")

    parser.add_argument("-alpha", type=float, default=90.0, help="Angle alpha (between b and c) in degree")
    parser.add_argument("-beta", type=float, default=90.0, help="Angle beta (between a and c) in degree")
    parser.add_argument("-gama", type=float, default=90.0, help="Angle gama (between a and b) in degree")

    args = parser.parse_args()

    a = args.a
    b = args.b
    c = args.c

    angl_bc = args.alpha * np.pi / 180.0
    angl_ac = args.beta * np.pi / 180.0
    angl_ab = args.gama * np.pi / 180.0

    L1, L2, L3 = FindLatVec(a, b, c, angl_bc, angl_ac, angl_ab)

    np.set_printoptions(precision=13)

    print( "L1 : ", L1 )
    print( "L2 : ", L2 )
    print( "L3 : ", L3 )


