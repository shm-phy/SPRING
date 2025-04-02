
import numpy as np


kB = 1.38064852E-23 #m2.kg.s-2.K-1

def NoseMass(N, T, tau) :

    Q = 3*(N-1)*kB*T*(tau**2)

    return Q

if __name__ == "__main__" :

    N = 250
    #T = 200 #K
    #T = 80 #K
    T = 300 #K
    
    dt = 2.0E-15 #fs
    Ndamp = 20
    tau = Ndamp*dt
    
    ERyd = 4.3597447222071E-18 / 2 #J
    fs2 = 1.0E-30 #s2
    MI = ERyd*fs2 #kg.m2
    
    Q = NoseMass(N, T, tau) #kg.m2
    
    print(Q/MI)
