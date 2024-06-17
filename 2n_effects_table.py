# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 15:37:16 2024

@author: aburtnerabt
"""
import pandas as pd
import numpy as np
import maplotlib.pyplot as plt

#TODO: Work through question 1, then use basis as utility class
data = pd.read_excel("C:/Users/aburtnerabt/Documents/K-State/Experimental Design/Assignment 2 Question 1 Data.xlsx")
data = data.sort_values(by=['Cooling Time','Glass Type','Molten-Glass Temperature'])
data

step1a = [data['Mold Temperature'].iloc[x]+data['Mold Temperature'].iloc[x+1] for
          x in range(0,int(data.shape[0]),2)]
step1b = [data['Mold Temperature'].iloc[x+1]-data['Mold Temperature'].iloc[x] for
          x in range(0,int(data.shape[0]),2)]

step1 = step1a.copy()
step1.extend(step1b)

def yates_step(data:list):
    """
    creates a list of equal length input list by first addign successive pairs
    and then subtracting sucessive pairs
    """
    #generate successive pair sums
    sums = [data[x]+data[x+1] for x in range(0,len(data),2)]
    
    #generate successive differences
    diffs = [data[x+1]-data[x] for x in range(0,len(data),2)]
    
    #extend sums
    sums.extend(diffs)
    return sums

def generate_standard_table(n, columns):
    """
    Generates standard format table for 2^n factorial experiment laid
    out in standard form
    """
    #get length of data frame
    length = 2**n
    
    #for every power of two, create binary arrays
    arrays = []
    for i in range(n):
        
        #instantiate empty factor list
        factor_levels = []
        
        #while factor level less than 2**n
        while len(factor_levels) < length:
            
            #append 2**i 0s
            for j in range(2**i):
                factor_levels.append(0)
            
            #append 2**i 1s
            for j in range(2**i):
                factor_levels.append(1)
        #add factor levels to arrays
        arrays.append(factor_levels)
    arrays = np.array(arrays).T
    table = pd.DataFrame(data=arrays, columns=columns)
    return table
"""
generate_standard_table(3, ['A','B','C'])

step_1 = yates_step(data['Mold Temperature'].tolist())
step_2 = yates_step(step_1)
step_3 = yates_step(step_2)
[x/8 if x == step_3[0] else x/4 for x in step_3]
[x / 8**.5 for x in step_3]


data = [71, 61, 90, 82, 68, 61, 87, 80, 61, 50, 89, 83, 59, 51, 85, 78]
step1 = yates_step(data)
step2 = yates_step(step1)
step3 = yates_step(step2)
step4 = yates_step(step3)
result = [x / 4 for x in step4]
"""
if __name__ == "__main__":
    pass
