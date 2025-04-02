
import numpy as np

np.set_printoptions(precision=13)

def CalculateMassVariance( Mass_arr, Abundance_frac ) :

    if ( np.sum(Abundance_frac) != 1.0 ) :
        print( "WARNING: Abundance_frac does not sum up to 1: ", np.sum(Abundance_frac) )

    AverageMass = np.sum( Mass_arr * Abundance_frac )

    print( "AverageMass: ", AverageMass )

    g2 = np.dot(Abundance_frac, (Mass_arr-AverageMass)**2 ) / (AverageMass**2)

    return g2

#Mass_arr = np.array([38.9637064848, 40.9618252561, 39.963998165], dtype=np.double) #K
#Abundance_frac = np.array([0.932581, 0.067302, 1.17E-4], dtype=np.double) #K

#Mass_arr = np.array([179.94747, 180.94800], dtype=np.double) #Ta
#Abundance_frac = np.array([0.0001176, 0.9998824], dtype=np.double) #Ta

#Mass_arr = np.array([69.9242474, 71.9220758, 72.9234589, 73.9211778, 75.9214026], dtype=np.double) #Ge
#Abundance_frac = np.array([0.2038, 0.2731, 0.0776, 0.3672, 0.0783], dtype=np.double) #Ge

#Mass_arr = np.array([73.9224764, 75.9192136, 76.9199140, 77.9173091, 79.9165213, 81.9166994], dtype=np.double) #Se
#Abundance_frac = np.array([0.0089, 0.0937, 0.0763, 0.2377, 0.4961, 0.0873], dtype=np.double) #Se

#Mass_arr = np.array([27.97692653442, 28.97649466434, 29.973770137], dtype=np.double) #Si
#Abundance_frac = np.array([0.92223, 0.04685, 0.03092], dtype=np.double) #Si

#Mass_arr = np.array([106.9050915, 108.9047558], dtype=np.double) #Ag
#Abundance_frac = np.array([0.518, 0.482], dtype=np.double) #Ag

#Mass_arr = np.array([203.9730436, 205.9744653, 206.9758969, 207.976652], dtype=np.double) #Pb
#Abundance_frac = np.array([0.0140, 0.241, 0.221, 0.524], dtype=np.double) #Pb

#Mass_arr = np.array([119.90402, 121.9030439, 122.9042700, 123.9028179, 124.9044307, 125.9033117, 127.9044631, 129.9062244], dtype=np.double) #Te
#Abundance_frac = np.array([0.0009, 0.0255, 0.0089, 0.0474, 0.0707, 0.188, 0.317, 0.341], dtype=np.double) #Te

#Mass_arr = np.array([12.0, 13.003354835336], dtype=np.double) #C
#Abundance_frac = np.array([0.9884, 0.0096], dtype=np.double) #C

#Mass_arr = np.array([27.9769265350, 28.9764946653, 29.973770137], dtype=np.double) #Si
#Abundance_frac = np.array([0.92223, 0.04685, 0.03092], dtype=np.double) #Si

#Mass_arr = np.array([10.01294, 11.00931], dtype=np.double) #B
#Abundance_frac = np.array([0.1965, 0.8035], dtype=np.double) #B

#Mass_arr = np.array([14.003074, 15.000109], dtype=np.double) #N
#Abundance_frac = np.array([0.996, 0.004], dtype=np.double) #N

#Mass_arr = np.array([34.9688527, 36.9659026], dtype=np.double) #Cl
#Abundance_frac = np.array([0.7576, 0.2424], dtype=np.double) #Cl

#Mass_arr = np.array([63.9291422, 65.9260334, 66.9271273, 67.9248442, 69.9253193], dtype=np.double) #Zn
#Abundance_perc = np.array([49.2, 27.7, 4, 18.5, 0.6], dtype=np.double) #Zn
#Abundance_frac = Abundance_perc / 100.0 #Zn

#Mass_arr = np.array([31.9720711744, 32.9714589099, 33.96786701, 35.96708070], dtype=np.double) #S
#Abundance_perc = np.array([94.8, 0.760, 4.37, 0.02], dtype=np.double) #S
#Abundance_frac = Abundance_perc / 100.0 #S

#Mass_arr = np.array([14.003074004251, 15.000108898266], dtype=np.double ) #N
#Abundance_frac = np.array([0.99578, 0.00337], dtype=np.double) #N

Mass_arr = np.array([15.994914619257, 16.999131755953, 17.999159612136], dtype=np.double ) #0
Abundance_frac = np.array([0.99738, 0.000367, 0.00187], dtype=np.double) #0

#Mass_arr = np.array([1.007825031898, 2.014101777844], dtype=np.double ) #H
#Abundance_frac = np.array([0.99972, 0.00028], dtype=np.double) #H

g2 = CalculateMassVariance( Mass_arr, Abundance_frac )

print( "g2 = ", g2)

