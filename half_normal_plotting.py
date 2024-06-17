# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 09:43:41 2024

@author: aburtnerabt
"""

import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.stats import norm
from scipy.special import ndtri
import pandas as pd
import math

#TODO create class for half normal plot
class HalfNormPlot:
    """
    Class that generates a half normal plot and provides functionality to
    estimate sigma hat using three different methods as outlined in lesson 2.
    
    Inputs
        data: pandas dataframe containing data column and label column
        data_col: str identifying which column is the data column
        label_col: str identifying which column is the label column
    """
    def __init__(self, 
                 data,
                 data_col,
                 label_col):
        
        
        self.half_norm_data = data
        self.data_col = data_col
        self.label_col = label_col
        self.abs_col = f"abs_{data_col}"
        
        #get abs value and sort
        self.half_norm_data[self.abs_col] = self.half_norm_data[data_col].apply(lambda x: abs(x))
        self.half_norm_data = self.half_norm_data.sort_values(by=self.abs_col)
        self.half_norm_data.reset_index(drop=True, inplace=True)
            
        #calculate Ri
        self.half_norm_data['r_i'] = self.half_norm_data[self.abs_col].rank(method="first")
        
        #calcualte Ri*
        self.half_norm_data['r_i*'] = (self.half_norm_data['r_i'] - .5) /data.shape[0]
        
        #calculate Pi
        self.half_norm_data['p_i'] = (self.half_norm_data['r_i*'] + 1) / 2
        
        #calculate Vi using
        self.half_norm_data['v_i'] = ndtri(self.half_norm_data['p_i'])
        
        return
    
    def plot_half_norm(self, 
                       title:str, 
                       data_percent=.95,
                       num_adjust=0):
        """
        Plots V_i against |X_i| to create the half normal plot
        
        Inputs
            title: str, text to place at top of title
            data_percent: float, percent of data to include in creation of
                                reference line
            num_adjust: number of data points to remove to account for 
                        significant effects
        """
        #create figure and axis
        fig, ax = plt.subplots(1,1, figsize=(8,8))
        
        #scatter plot the data points
        ax.scatter(self.half_norm_data[self.abs_col], 
                   self.half_norm_data['v_i'])
        ax.set_xlabel(self.abs_col)
        ax.set_ylabel("V_i")
        ax.set_title(title)
        
        #add line through first 90% of data points for rough identification
        #of outliers
        
        #if num adjust is true, we use all the remaining data, otherwise use
        #a (large) fraction of the data to fit the regression
        #also instantiate reference line data with a point at 0 for better
        #plotting of half normal plots
        self.ref_line_data = [0]
        self.ref_line_data.extend(self.half_norm_data[self.abs_col].tolist())
        if num_adjust == 0:
            
            #get number of data points to use to fit line
            n = math.floor(data_percent*self.half_norm_data.shape[0])
            
            #make X and Y arrays
            X = self.half_norm_data[self.abs_col].iloc[:n+1]
            y = self.half_norm_data['v_i'].iloc[:n+1]
            
            #create model and fit
            self.model = sm.OLS(y, X, hasconst=False)
            self.model_results = self.model.fit()

            #calcualte and plot reference line
            ref_line = [x*self.model_results.params[self.abs_col]  for x 
                        in self.ref_line_data]
            ax.plot(self.ref_line_data, ref_line, color='red')
            
        else:
            #get number of data points to fit
            n = self.half_norm_data.shape[0]-num_adjust
            
            #make X and Y arrays
            self.X = self.half_norm_data[self.abs_col].iloc[:n]
            y = self.half_norm_data['v_i'].iloc[:n]
            
            #create model and fit model
            self.model_adj = sm.OLS(y, self.X, hasconst=False)
            self.model_adj_results = self.model_adj.fit()
            
            #plot reference line
            ref_line = [x*self.model_adj_results.params[self.abs_col]  for x 
                        in self.ref_line_data]
            ax.plot(self.ref_line_data, ref_line, color='red')
            
            #add labels to plot
            self.label_points = [(x, y) for x, y in 
                            zip(self.half_norm_data[self.abs_col].iloc[-num_adjust:],
                                self.half_norm_data['v_i'].iloc[-num_adjust:])]
            self.labels = self.half_norm_data[self.label_col].iloc[-num_adjust:]
            
            for label, point in zip(self.labels, self.label_points):
                text_loc = (point[0]*.95, point[1]*1.05)
                ax.annotate(label, point, xytext=text_loc, size=12)
        
        return ax
    
    def sigma_closest_to_one(self):
        """
        Function that returns rough estimate of sigma using the data point
        that is closest to one
        """
        #get distance from v_i and 1
        dist_from_one = [abs(x-1) for x in self.half_norm_data['v_i']]
        
        #get the minimum of distance to one
        min_dist = min(dist_from_one)
        loc = dist_from_one.index(min_dist)
        return self.half_norm_data[self.abs_col].iloc[loc]
    
    def sigma_from_regression(self):
        """
        Function that gives estimate of signma using the method of
        dividing 1 by the slope of the fitted regression line 
        """
        return 1/self.model_results.params[self.abs_col]
    
    def sigma_from_adjusted_regression(self):
        return 1/self.model_adj_results.params[self.abs_col]

        

"""
#read in table 4.1 from "Analysis of Messy Data Volume 2"
data = pd.read_csv('C:/Users/aburtnerabt/Documents/K-State/Experimental Design/table_4.1.csv')

#create half normal plot
table_4_1_test = HalfNormPlot(data, data_col='X', label_col='Label')

#compare to outcome in lecture slides or table 4.1
table_4_1_test.half_norm_data

#compare to Fig 4.1
table_4_1_test.plot_half_norm("Half normal plot of data in table 4.1", 
                              data_percent=.95)

#compare to sigma hat estimate page 67
table_4_1_test.sigma_closest_to_one()

#compare to sigma hat estimate lecture slides example 3 or page 68
table_4_1_test.sigma_from_regression()

#below executions have no validtion data in book or lecture
table_4_1_test.plot_half_norm("Adjusted half norm plot", num_adjust=1)
table_4_1_test.sigma_from_adjusted_regression()
table_4_1_test.half_norm_data.shape


test = np.random.rand(20,1)
(test - .5) /10
test_sort = np.sort(test,axis=0)
type(test) == np.ndarray

np.array(range(1,test.shape[0]+1))


half_norm = HalfNormPlot(test)
half_norm = HalfNormPlot(test.tolist())
norm.pdf(test)
"""

if __name__ == "__main__":
    pass
