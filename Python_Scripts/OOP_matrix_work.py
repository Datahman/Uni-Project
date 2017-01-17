# Switch IX Hamiltonian to a Class type.

import numpy as np
from scipy import linalg as LA
from math import pi
import cmath
from datetime import datetime
from scipy.constants import hbar, m_e, eV, c
from scipy.optimize import fsolve
from scipy import optimize
import pandas as pd
import matplotlib.pyplot as plt
class Matrix(object):
    
    def __init__(self,m,n,BMagnetic,deltab,k,thetak):
        self.m = m
        self.n = n
        self.BMagnetic = BMagnetic
        self.deltab=deltab
        self.k = k
        self.thetak = thetak
        self.m_e = m_e
        # For the moment have all the Hamiltonian constants accessible to class-wide structure.
        self.mub = mub = 57.88
        self.mm_e =  (self.m_e * c**2) /(1.6*(10**-25)) # Electron mass.
        self.hbarc = (1.24 * (10**6)) / (2 * pi) 
        self.m_ex = (0.21) * (self.mm_e)    
        self.gemub = -0.01 * self.mub
        self.ghmub = 0.0085 * self.mub
        self.alphae = 0
        self.betae = 2.7
        self.alphah = 0
        self.betah = 0.92
        self.gammah = 0.0
        self.EbmEd = 5
        self.deltad = -13
       
    def make_hhh1(self): # This method creates the on-site diagonal matrix elements, and phase associated with the momentum.
         
        self.phasef = cmath.cos(self.thetak)+ 1j * cmath.sin(self.thetak) 
        self.phasefc = cmath.cos(self.thetak) - (1j *  cmath.sin(self.thetak))
        self.hhh1 = float (((self.hbarc * self.k)**2) / (2*self.m_ex))
        return self.hhh1            # Output: 72.4909719004
                                                                #IX mass: m_ex= 1.07455761652e+11 , hbarc = 197352.129434, k =20
    
    def Make_matrix(self):
        self.mylist = []
        
        for column in range(1,self.m+1):
            for row in range(1,self.n+1):
                if column==row: # For diagonal elements
                    if column==1:
                        #self.mylist.append(1)
                        self.mylist.append(((self.gemub * self.BMagnetic)  / 2)+ self.make_hhh1()) 
                    elif column==2:
                        self.mylist.append(((self.gemub * self.BMagnetic)  / 2)+ self.EbmEd/2)  
                        
                    elif column==3:
                        #self.mylist.append(24)
                        self.mylist.append(((-self.gemub * self.BMagnetic)  / 2)- self.EbmEd/2)  
                    elif column==4:
                        #self.mylist.append(0)
                        self.mylist.append((( self.gemub * self.BMagnetic )/2) - self.EbmEd/2)   
                       
                      
                else:
                    if column != row:
                        if column ==1 and row ==2: 
                            #self.mylist.append(0)
                            self.mylist.append(-self.deltab)
                        elif column ==1 and row ==3:
                            #self.mylist.append(0)
                            self.mylist.append((self.alphae * self.phasefc * self.k) + (self.betae * self.phasef * self.k ))
                        elif column ==1 and row ==4:
                            #self.mylist.append(0)
                            self.mylist.append((self.alphah * self.phasefc * self.k) + (self.betah * self.phasef * self.k  ) + ((self.gammah)* (self.phasefc * self.k)**3))
                            
                        elif column ==2 and row ==1: 
                            #self.mylist.append(0)
                            self.mylist.append(-self.deltab)
                        elif column ==2 and row ==3:
                            #self.mylist.append(0)
                            self.mylist.append((self.alphah * self.phasef * self.k ) + (self.betah * self.phasefc * self.k) + (self.gammah * (self.phasef * self.k)**3))
                        elif column ==2 and row ==4:
                            #self.mylist.append(0)
                            self.mylist.append((self.alphae * self.phasef * self.k ) + (self.betae  * self.phasefc * self.k))
                            
                        elif column ==3 and row ==1: 
                            #self.mylist.append(3)
                            self.mylist.append((self.alphae * self.phasef * self.k ) + (self.betae * self.phasefc * self.k ))
                        elif column ==3 and row ==2:
                            #self.mylist.append(0)
                            self.mylist.append((self.alphah * self.phasefc * self.k ) + (self.betah * self.phasef * self.k) + (self.gammah * (self.phasefc*self.k)**3))
                        elif column ==3 and row ==4:
                            #self.mylist.append(3)
                            self.mylist.append(-self.deltad)
                        
                        elif column ==4 and row ==1: 
                            #self.mylist.append(0)
                            self.mylist.append((self.alphah * self.phasef * self.k) + (self.betah * self.phasefc * self.k ) +(self.gammah* (self.phasef * self.k)**3))
                        elif column ==4 and row ==2:
                            #self.mylist.append(23)
                            self.mylist.append((self.alphae * self.phasefc * self.k) + (self.betae * self.phasef * self.k))
                        elif column ==4 and row ==3:
                            #self.mylist.append(0)
                            self.mylist.append(-self.deltad)
                            
        matrix = np.array([self.mylist],dtype='complex64')
        matrix = matrix.reshape(4,4)
        return matrix   
        
    def Diag_matrix(self): #Self-note: Member function calling on instance self.
        self.eigval_list=[]
        self.eigvec_list = []
        eigvals,eigvecs = LA.eig(self.Make_matrix())
        self.eigval_list.append(eigvals)
        self.eigvec_list.append(eigvecs)
        return self.eigval_list, self.eigvec_list #eigvals,eigvecs
    

        
class Fermi_Momenta(Matrix):
    def __init__(self,m,n,BMagnetic,deltab,k,thetak,NumkGridPoints,k_max,e_F):
        Matrix.__init__(self,m,n,BMagnetic,deltab,k,thetak)
        self.NumkGridPoints = NumkGridPoints
        self.k_max = k_max
        self.e_F = e_F
    
    def Diag_matrix(self):
        a = Matrix.Diag_matrix(self) # Personal note: If no argument provided then TypeError will be thrown, to fix use Matrix class instance "self". 
        b = list(a)
        self.k_list = []
        self.eigval_list = []
        for k_index in range(1,self.NumkGridPoints+1):
            k_vals = (self.k_max * (1-k_index)) / (1- self.NumkGridPoints) # Grid point setup for momentum(k) values.
            self.k_list.append(float(k_vals))
        #print b
        self.eigval_list.append(b[0])
             
        return self.eigval_list, self.k_list
    def All_k(self): # Method to get eigenvalues for all k's as retreived from Diag_matrix (subclass)
        self.Diag_matrix()
        # Each state eigen values at every k.
        self.eigN1 = []   
        self.eigN2 = []
        self.eigN3 = []
        self.eigN4 = []
        for each_k in self.k_list:
            eigval,eigvec = Matrix(4,4,0.5,0.5,each_k,10).Diag_matrix() 
            for each_state_val in eigval:
                 self.eigN1.append(each_state_val[0]),self.eigN2.append(each_state_val[1])
                 self.eigN3.append(each_state_val[2]),self.eigN4.append(each_state_val[3])   
                 #print each_state_val[0],each_state_val[1],each_state_val[2],each_state_val[3]
                 
        return self.eigN1,self.eigN2,self.eigN3,self.eigN4       
        #super(Fermi_Momenta,self).Diag_matrix()
    
    def Dispersion_Plot(self):
        self.Diag_matrix() #  Personal note: Call to subclass method to use self.k_list.
        self.All_k()
        
        plt.plot(self.k_list,(self.eigN1),'g-', label='State 1', linewidth=2)
        plt.plot(self.k_list,(self.eigN2),'b-', label='State 2', linewidth=2)
        plt.plot(self.k_list,(self.eigN3),'r-', label='State 3', linewidth=2)
        plt.plot(self.k_list,(self.eigN4),'y-', label='State 4', linewidth=2)
        #plt.ylabel('est')
        plt.title('Dummy graph')
        plt.ylabel('Energy')
        plt.xlabel('kF')
        plt.legend()
        plt.show()
    

a = Matrix(4,4,0.5,0.5,20,10).Make_matrix()
b = Matrix(4,4,0.5,0.5,20,10).Diag_matrix()
c_ = Matrix(4,4,0.5,0.5,20,10)
d = c_.make_hhh1()
e = Fermi_Momenta(4,4,0.5,0.5,20,10,10,20,20).Diag_matrix() #In the following order: {m,n,BMagnetic,deltab,k,thetak,NumkGridPoints,k_max,e_F}
f = Fermi_Momenta(4,4,0.5,0.5,20,10,100,50,40).Dispersion_Plot()
print "#######",f

#for states in range(0,self.m): # To do.. make the diag matrix return eigen values. Done !
# To do: Add hhh1 method return values to subclass Fermi_Momenta. Done !
# To do: Check the eigval-momentum graph. Done.
# To do: Fix eigenvalue graph. 
# To do: Maybe later change subclass Diagmatrix name.
 