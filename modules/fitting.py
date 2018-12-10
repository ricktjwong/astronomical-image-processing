# -*- coding: utf-8 -*-
"""
Created on Mon Dec 03 11:19:30 2018

@author: drs16
"""

import numpy as np
import matplotlib.pyplot as plt

def sersic(Ie, Re, R, n):
    """ Ie is intensity at half light radius Re; Ie, Re are constant parameters. """
    b = 2*n**(-1/3)
    return Ie * np.exp(-b * ( (R/Re)**(1/n) - 1) )
    

class Sersic_fitter:
    def __init__ (self, Ie, Re, data_R, data_I):
        self.data_R = data_R
        self.data_I = data_I
        self.Ie = Ie
        self.Re = Re
        self.model_I = None
        self.model_continuous = None
        self.R_ = None
        self.res_sum_sq, self.reg_sum_sq, self.tot_sum_sq = 0,0,0
        self.chi_sq = 0
        self.R_sq = 0
        self.best_n = None
      
    def get_n(self):
        return self.best_n
        
    def get_cont(self):
        return self.model_continuous
  
    def sersic(self, R, n):
            b = 2*n**(-1/3)
            return self.Ie * np.exp(-b * ( (float(R)/float(self.Re)) ** (1/n) - 1) )

    def function_fit(self,n):
        self.model_I = [self.sersic(i,n) for i in self.data_R]
    
    def function_cont(self,n, numOfPoints):
        R_ = np.linspace(self.data_R[0], self.data_R[-1], numOfPoints)
        self.R_ = R_
        self.model_continuous = [self.sersic(i,n) for i in R_]
    
    def res_sum_square(self):
        for i,j in zip(self.model_I, self.data_I):
            residual_sq = (i - j) * (i - j)
            self.res_sum_sq += residual_sq  
                
    def chi_square(self):
        chi_sum = 0
        for i,j in zip(self.model_I, self.data_I):
            current_chi = ((i - j) * (i - j))/i
            chi_sum += current_chi
        self.chi_sq = chi_sum
        
    def total_sum_square(self):
        mean = np.mean(self.data_I)
        cumulative_sum = 0
        for i in self.data_I:
            cumulative_sum += (i - mean) * (i - mean)
        self.tot_sum_sq = cumulative_sum
        
    def reg_sum_square(self):
        mean = np.mean(self.data_I)
        cumulative_sum = 0
        for i in self.model_I:
            cumulative_sum += (i - mean) * (i - mean)
        self.reg_sum_sq = cumulative_sum
        
    def R_square(self):
        self.function_fit(0.75)
        self.reg_sum_square()
        self.total_sum_square()
        R = (self.reg_sum_sq/self.tot_sum_sq)
        self.R_sq = R
    
    def sweep_optimise(self):
        n_list = []
        sum_sq_list = []
        for i in np.arange(0.5, 10, 0.05):
            self.res_sum_sq = 0                 # reset to zero each iteration
            self.function_fit(i)
            self.res_sum_square()
            n_list.append(i)
            sum_sq_list.append(self.res_sum_sq)
        self.best_n = n_list[np.argmin(sum_sq_list)]
        

            
            
R = np.arange(3,13,1)
#I = [46278, 63769, 72319, 46935, 37199, 20208, 13447, 11850, 8649, 6718, 4348, 3514]
I = [72319, 46935, 37199, 20208, 13447, 11850, 8649, 6718, 4348, 3514]
Ie = 8649
Re = 9.194
fitter = Sersic_fitter(Ie, int(Re), R, I)
fitter.sweep_optimise()
n = fitter.best_n
fitter.function_cont(0.75,50)
fitter.function_fit(0.75)
fitter.R_square()

print (fitter.R_sq)
print (n)
""" Plotting """
plt.figure()
plt.plot(fitter.R_, fitter.model_continuous)
plt.plot(R, I, 'r.')
plt.show()

