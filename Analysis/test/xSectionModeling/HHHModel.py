from sympy import *
from numpy import matrix
from numpy import linalg
from sympy import Matrix

class HHHSample:
    def __init__(self, val_k3, val_k4, val_xs, label):
        self.val_k3  = val_k3
        self.val_k4  = val_k4
        self.val_xs  = val_xs
        self.label   = label

####################
class HHHFormula:
    def __init__(self, sample_list):
        self.sample_list = sample_list
        self.build_matrix()
        self.sigmaEval=None
        self.calculatecoeffients()

    def build_matrix(self):
        """ create the matrix M in this object """

        if len(self.sample_list) != 9:
            print ("[ERROR] : expecting 9 samples in input")
            raise RuntimeError("malformed vbf input sample list")
        M_tofill = [
            [None,None,None,None,None,None,None,None,None],
            [None,None,None,None,None,None,None,None,None],
            [None,None,None,None,None,None,None,None,None],
            [None,None,None,None,None,None,None,None,None],
            [None,None,None,None,None,None,None,None,None],
            [None,None,None,None,None,None,None,None,None],
            [None,None,None,None,None,None,None,None,None],
            [None,None,None,None,None,None,None,None,None],
            [None,None,None,None,None,None,None,None,None]
        ]
        
        self.xSections={}
        for isample, sample in enumerate(self.sample_list):
            print(isample, " k3, k4  = ", sample.val_k3, sample.val_k4, sample.val_xs)

            ## implement the 9 scalings
            M_tofill[isample][0] = sample.val_k3**4                             
            M_tofill[isample][1] = sample.val_k3**3
            M_tofill[isample][2] = sample.val_k3**2
            M_tofill[isample][3] = sample.val_k3**2 * sample.val_k4
            M_tofill[isample][4] = sample.val_k3    
            M_tofill[isample][5] = sample.val_k3 * sample.val_k4
            M_tofill[isample][6] = sample.val_k4**2
            M_tofill[isample][7] = sample.val_k4
            M_tofill[isample][8] = 1.0

            self.xSections['s'+str(isample+1)]=sample.val_xs


        #print (M_tofill)
        #print("Input Matrix defined as  :\n")
        #for row in  M_tofill:
        #    print(end="   ")
        #    for item in row:
        #        print("{0:0.3f}".format(item), " | ",end="")
        #    print()
 
        self.M = Matrix(M_tofill)
        #print ("\n Determinant of the matrix : " , det(self.M) )

    def calculatecoeffients(self):
        """ create the function sigma and the nine coefficients in this object """

        try: self.M
        except AttributeError: self.build_matrix()
        ##############################################    
        k3, k4, A, B, C, D, E, F, G, H, I, s1, s2, s3, s4, s5, s6, s7, s8, s9 = symbols('k3, k4, A, B, C, D, E, F, G, H, I, s1, s2, s3, s4, s5, s6, s7, s8, s9')
        ### the vector of couplings
        ### the vector of couplings
        c = Matrix([
            [k3**4 ] ,
            [k3**3] ,
            [k3**2] ,
            [k3**2 * k4] ,
            [k3] ,
            [k3 * k4] ,
            [k4**2] ,
            [k4] ,
            [1] ,
        ])
        ### the vector of components
        v = Matrix([
            [A] ,
            [B] ,
            [C] ,
            [D] ,
            [E] ,
            [F] ,
            [G] ,
            [H] ,
            [I] 
        ])
        ### the vector of samples (i.e. cross sections)
        s = Matrix([
            [s1] ,
            [s2] ,
            [s3] ,
            [s4] ,
            [s5] ,
            [s6] ,
            [s7] ,
            [s8] ,
            [s9] 
        ])
        ####    
        Minv   = self.M.inv()
        self.coeffs = c.transpose() * Minv # coeffs * s is the sigma, accessing per component gives each sample scaling
        self.sigma  = self.coeffs*s
        substitutions=[]
        for ky in self.xSections:
            substitutions.append((ky,self.xSections[ky]))
        self.sigmaEval=self.sigma.subs(substitutions)

    def evaluateSigma(self,params):
        substitutions=[]
        if self.sigmaEval==None:
            print("Evaliate the matrix first !! ")
            return -1.0
        for ky in params:
            substitutions.append((ky,params[ky]))
        return self.sigmaEval.subs(substitutions)[0]


####################
# set of all cross sections

txt="""
1  |  k3 =  0.000  	 k4 =  0.000 	 xSec =  9.605
2  |  k3 =  0.000  	 k4 =  1.000 	 xSec =  9.967
3  |  k3 =  1.000  	 k4 =  0.000 	 xSec =  3.599
4  |  k3 =  1.000  	 k4 =  1.000 	 xSec =  3.253
5  |  k3 =  2.000  	 k4 =  1.000 	 xSec =  2.539
6  |  k3 =  2.000  	 k4 =  3.000 	 xSec =  1.403
7  |  k3 =  3.000  	 k4 =  0.000 	 xSec =  5.059
8  |  k3 =  3.000  	 k4 =  2.000 	 xSec =  3.570
9  |  k3 =  -0.500  	 k4 =  0.500 	 xSec =  17.120
"""
print("xSections Used \n",txt)
HHH_List=[]
HHH_List.append(HHHSample(0.0 , 0.0 ,val_xs=9.605, label='HHH_k3_0_k4_0'))
HHH_List.append(HHHSample(0.0 , 1.0 ,val_xs=9.967, label='HHH_k3_0_k4_1'))
HHH_List.append(HHHSample(1.0 , 0.0 ,val_xs=3.599, label='HHH_k3_1_k4_0'))
HHH_List.append(HHHSample(1.0 , 1.0 ,val_xs=3.253, label='HHH_k3_1_k4_1'))
HHH_List.append(HHHSample(2.0 , 1.0 ,val_xs=2.539, label='HHH_k3_2_k4_1'))
HHH_List.append(HHHSample(2.0 , 3.0 ,val_xs=1.403, label='HHH_k3_2_k4_3'))
HHH_List.append(HHHSample(3.0 , 0.0 ,val_xs=5.059, label='HHH_k3_3_k4_0'))
HHH_List.append(HHHSample(3.0 , 2.0 ,val_xs=3.570, label='HHH_k3_3_k4_2'))
HHH_List.append(HHHSample(-0.50 , 0.5 ,val_xs=17.1200,label='HHH_k3_m0p5_k4_0p5'))

print("Making the Object with ",len(HHH_List)," samples \n") 
gHere = HHHFormula(HHH_List)
print("Formula =  ",gHere.sigma[0])


