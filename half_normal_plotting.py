# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 09:43:41 2024

@author: aburtnerabt
"""

import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.stats import norm
import pandas as pd

#TODO create class for half normal plot
class HalfNormPlot:
    """
    Class that generates a half normal plot and provides functionality to
    estimate sigma hat using three different methods as outlined in lesson 2.
    
    Inputs
        data: np.array, must be 2D array of shape [n,1]
    """
    def __init__(self, data):
        
        #validate to be (2,1) numpy array
        assert type(data) == np.ndarray, f"Expected data to be numpy array, got {type(data)}"
        assert data.shape[1] == 1, f"Expected single column of data, got passed {data.shape[1]} columns"
        
        self.half_norm_data = pd.DataFrame()
        
        #get abs value and sort
        self.half_norm_data['abs_x'] = np.sort(np.absolute(data),axis=0)
            
        #calculate Ri
        self.half_norm_data['r_i'] = np.array(range(1,data.shape[0]+1))
        
        #calcualte Ri*
        self.half_norm_data['r_i*'] = (self.half_norm_data['r_i'] - .5) /10
        
        #calculate Pi
        self.half_norm_data['p_i'] = (self.half_norm_data['r_i*'] + 1) / 2
        
        #calculate Vi
        self.half_norm_data['v_i'] = norm.pdf(self.half_norm_data['p_i'])
        

            



#TODO calculate Ri, Ri*, Pi, Vi



test = np.random.rand(20,1)
(test - .5) /10
test_sort = np.sort(test,axis=0)
type(test) == np.ndarray

np.array(range(1,test.shape[0]+1))


half_norm = HalfNormPlot(test)
half_norm = HalfNormPlot(test.tolist())
norm.pdf(test)
