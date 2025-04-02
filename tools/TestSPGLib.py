
import spglib
import numpy as np

# Mind that the a, b, c axes are given in row vectors here,
# but the formulation above is given for the column vectors.

PREC = 1.0E-6

#lattice = [[4.0659928322, 0.0000000000, 0.0000000000],  # a
#           [2.0329964161, 3.5212530843, 0.0000000000],  # b
#           [2.0329964161, 1.1737510281, 3.3198692456]]  # c
#
##points = [[0.00, 0.00, 0.00],  #Ga
##          [0.25, 0.25, 0.25]]  #As
#points = [[0.25, 0.25, 0.25],  #Ga
#          [0.50, 0.50, 0.50]]  #As
#numbers = [1, 2]

#lattice = [[2.513241657,  0.000000000,  0.000000000],
#           [1.256620828,  2.176531121,  0.000000000],
#           [1.256620828,  0.725510374,  2.052053222]]
#points = [[0.50, 0.50, 0.50],  
#          [0.75, 0.75, 0.75]] 
#numbers = [1, 1]

lattice = [[3.85392288,     0.0000000000000, 0.0000000000000],
           [1.92696144,     3.3375951183061, 0.0000000000000],
           [1.92696144,     1.1125317061020, 3.1467148546791]]
points = [[0.25, 0.25, 0.25],  #Zn
          [0.50, 0.50, 0.50]]  #S
#points = [[0.00, 0.00, 0.00],  #Zn
#          [0.75, 0.75, 0.75]]  #S
numbers = [1, 2]


cell = (lattice, points, numbers)
dataset = spglib.get_symmetry_dataset(cell, hall_number=0)

print("=================================================")
print("Space group type: %s (%d)"
      % (dataset['international'], dataset['number']))
print("Transformation matrix:")
for x in dataset['transformation_matrix']:
    print("  %2d %2d %2d" % tuple(x))
print("Origin shift: %f %f %f" % tuple(dataset['origin_shift']))
print("=================================================")
print()

symmetry = spglib.get_symmetry(cell, symprec=PREC)
print("=================================================")
print("Number of Rotation Matrix :", len(symmetry['rotations']), "|| symprec=", PREC )
print("=================================================")
print()

num = 1
print("=================================================")
print("Printing", num, "th rotation and translation")
print("The rotation matrix considered :")
print(symmetry['rotations'][num])
print("Translation :")
print(symmetry['translations'][num])
print("=================================================")
print()

"""
R = symmetry['rotations'][num]
T = symmetry['translations'][num]

atm_typ = 1
atom = np.array(points[atm_typ])
print("Atom to operate on : (The atom index=", numbers[atm_typ], ")")
print(atom)
print("After the operation atom goes to :")
print(R.dot(atom)+T)

for i in range(len(symmetry['rotations'])) :
    print(symmetry['translations'][i])
    mat = symmetry['rotations'][i]
    if (mat == np.eye(mat.shape[0])).all() :
        print("Identity found :", i)
        print("Corresponding translation :", symmetry['translations'][i])

"""
