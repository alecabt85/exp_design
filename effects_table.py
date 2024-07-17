# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 15:37:16 2024

@author: aburtnerabt
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from collections import Counter
import string

def get_gen_interaction(effects):
    all_together = ''.join(effects)
    counter = Counter(all_together)
    residual = []
    for k,v in counter.items():
        if v % 2 == 1:
            residual.append(k)
    residual = sorted(residual)
    return ''.join(residual)

def get_exp_alias(factors, defining_effects):
    """
    Function to get alias sets for a fractional experiment by getting
    generalized interactions for all factors for defining effects and the 
    generalized interaction of the defining effect
    
    Inputs
    factors: list of strings of factors
    defining_effects: list of strings of defining effects
    """
    #get GI of Defining effects and add to list
    defining_effects = defining_effects.copy()
    def_eff_gi = get_gen_interaction(defining_effects)
    if len(defining_effects) >1:
        combos = []
        #iterate over all possible combinations and get gen interaction
        for i in range(1,len(defining_effects)+1):
            combos.extend([get_gen_interaction(x) for x in combinations(defining_effects,i)])
    
    defining_effects.extend(combos)
    
    defining_effects = list(set(defining_effects))
    
    #Build dataframe index as all combinations
    idx = []
    for i in range(1,len(factors)+1):
        idx.extend(combinations(factors,i))
    
    #make tuples strings
    idx = [''.join(x) for x in idx]
    idx.append('')
    
    #iterate over defining effects and factors and determin alias sets
    #make empty df to hold results
    alias_sets = pd.DataFrame(columns=defining_effects, index=idx)
    for col in defining_effects:
        for row in idx:
            alias_sets.loc[row,col] = get_gen_interaction([row,col])
            
    #make list of tuples
    alias_tuples = [tuple(sorted(x)) for x in alias_sets.itertuples()]
    
    #make unique list of tuples
    alias_tuples = list(set(alias_tuples))
    
    #sort tuples by length
    alias_tuples = [sorted(x, key=lambda x: len(x)) for x in alias_tuples]
    
    #make df
    alias_sets = pd.DataFrame(data=alias_tuples, columns=[f"alias_{i}" for i in range(1,len(alias_tuples[0])+1)])
    alias_sets = alias_sets.sort_values(by='alias_1')
    return alias_sets

def get_confounding_effects(defining_effects):
    """
    Generates the defining effects of an experiment by multiplying all 
    combinations of defining effects together
    """
    effects = {}
    for i in range(2,len(defining_effects)+1):
        pairs = combinations(defining_effects,i)
        for p in pairs:
            effects.update({'-'.join(p):get_gen_interaction(p)})
    return effects

def check_alias(int_ord, alias_sets, defining_effects):
    """
    Checks for aliasing of effects up to interaction order "int_ord" by 
    iterating over the "alias_sets" df and checking to see if any of the aliases
    are in the defining effects listing, thereby verifying if any effects
    are aliased with the defining effects
    """
    check = []
    for tup in alias_sets.itertuples():
        if len(tup.Index)>int_ord or len(tup.Index) < 1:
            continue
        check.extend(list(map(lambda x: x in defining_effects, tup)))
    
    return any(check)
        

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

def generate_standard_table(n):
    """
    Generates standard format table for 2^n factorial experiment laid
    out in standard form
    """
    #make columns
    columns = [x for x in string.ascii_uppercase[:n]]
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

def encode_effects(table, factors):
    #create empty list of columns
    columns = []
    
    
    #make +/-1 encoded version of table
    ctable = table.copy()
    
    for factor in factors:
        ctable[factor] = ctable[factor].map({0:-1,1:1})
    
    #make combinations of factors
    for i in range(1, len(factors)+1):
        columns.extend(list(combinations(factors,i)))
        
    #calculate contrast columns
    for col in columns:
        col_name = f"{''.join(col)} Effect"
        table[col_name] = ctable[[x for x in col]].prod(axis=1)
    return table
        

class Experiment2N:
    """
    Class to represent non-replicated 2^n experiments.  Intended work flow is
    as follows:
        1) Generate the standard format table
        2) Capture results in the standard format table and then feed 
            that table back into the "calculate_effects" method. Ensure that
            only main factor columns are in the results dataframe fed back into
            this method or index issues will happen
            
    For a half normal plot use the HalfNormPlot class
    """
    def __init__(self, 
                 factors,
                 block_size=0):
        self.n = len(factors)
        self.factors = factors
        self.block_size = block_size
        
        if block_size == 0:
            self.num_blocks = 1
        else:
            self.num_blocks = 2**(self.n) / self.block_size
        self.num_test_runs = 2**(self.n)
        
        return
    
    def generate_table(self, generate_contrast=True):
        """
        Generates standard format table for 2^n factorial experiment laid
        out in standard form
        """

        #for every power of two, create binary arrays
        arrays = []
        for i in range(self.n):
            
            #instantiate empty factor list
            factor_levels = []
            
            #while factor level less than 2**n
            while len(factor_levels) < self.num_test_runs:
                
                #append 2**i 0s
                for j in range(2**i):
                    factor_levels.append(0)
                
                #append 2**i 1s
                for j in range(2**i):
                    factor_levels.append(1)
            #add factor levels to arrays
            arrays.append(factor_levels)
        arrays = np.array(arrays).T
        self.table = pd.DataFrame(data=arrays, columns=self.factors)
        
        #if you don't want to see contrast columns return just the table
        if generate_contrast == False:
            return self.table.copy()
            
        
        #make table of -1/+1 to calculate contrasts
        self.contrast_table = self.table.copy()
        for col in self.contrast_table.columns:
            self.contrast_table[col] = self.contrast_table[col].map({0:-1, 1:1})
        
        #add contrast columns
        self.contrast_columns = []
        
        #generate all possible combos of factors
        for i in range(1, self.n+1):
            self.contrast_columns.extend(list(combinations(self.factors,i)))
            
        #for each combination calculate column and make it
        for col in self.contrast_columns:
            col_name = ''.join([x for x in col])
            self.table[f"{col_name}_contrast"] = self.contrast_table[[x for x in col]].prod(axis=1)
        return self.table.copy()
    
    def calculate_effects(self, results: pd.DataFrame, response_col: str):
        """
        Function that calculates the main and interaction effects
        of the 2n experiment.  This does NOT calculate block effects and 
        assumes there is one complete replication.  
        
        The data frame MUST HAVE the factors in the same column order as the
        class experiment table, the columns should be 0/1 encoded for high-low
        treatment, and the table should be in standard order.
        
        The ONLY columns that should be in the results df are the factors,
        any contrast columns should be dropped out.  The results df gets
        joined to the existing experiment table using factors as join keys,
        any column outside of the factors such as "a_contrast" will cause
        indexing issues.
        """
        #set index of results to the factor columns
        results = results.set_index(self.factors)
        
        #create results table by merge on factor column index
        self.effects_table = self.table.copy().set_index(self.factors)
        self.effects_table = self.effects_table.merge(results, 
                                                       how='left', 
                                                       left_index=True, 
                                                       right_index=True)
        
        #reset the index
        self.effects_table = self.effects_table.reset_index()
        
        #make dataframe holder for calculations
        self.calculated_effects = pd.DataFrame(data=np.zeros((2,len(self.contrast_columns))),
                                               columns=[f"{''.join(x)}_effect" for x in self.contrast_columns],
                                               index=['effect','std_effect'])
        
        #calculate each interaction effect and make
        for col in self.contrast_columns:
            col_name = ''.join(col)
            con_col = f"{col_name}_contrast"
            effect_col = f"{col_name}_effect"
            effect = sum([x*y for x,y in zip(self.effects_table[con_col],
                                             results[response_col])])
            self.calculated_effects.at['effect',effect_col] = effect
            self.calculated_effects.at['std_effect',effect_col] = effect / (len(
                self.effects_table[con_col])**.5)
            
        return self.calculated_effects.copy()
    
def generate_standard_table3N(n):

    #make columns
    columns = [x for x in string.ascii_uppercase[:n]]
    #get length of data frame
    length = 3**n
    
    #for every power of two, create binary arrays
    arrays = []
    for i in range(n):
        
        #instantiate empty factor list
        factor_levels = []
        
        #while factor level less than 3**n
        while len(factor_levels) < length:
            
            #append 2**i 0s
            for j in range(3**i):
                factor_levels.append(0)
            
            #append 2**i 1s
            for j in range(3**i):
                factor_levels.append(1)
                
            for k in range(3**i):
                factor_levels.append(2)
        #add factor levels to arrays
        arrays.append(factor_levels)
    arrays = np.array(arrays).T
    table = pd.DataFrame(data=arrays, columns=columns)
    
    #encode contrast columns
    linear_c = {
        0:-1,
        1: 0,
        2: 1}
    quad_c = {
        0: 1,
        1: -2,
        2: 1}
    
    #for each column create the linear and quadratic contrast
    for col in columns:
        table[f"{col}_L"] = table[col].map(linear_c)
        table[f"{col}_Q"] = table[col].map(quad_c)
        
    #get all contrast columns
    contast_cols = [x for x in table.columns if "_" in x]
    
    #get all linear and quadratic columns and make all combinations
    combos = []
    for i in range(2, n+1):
        combos.extend(combinations(contast_cols,i))
    combo_columns = [" x ".join(x) for x in combos]
    
    #zip together combo columns and combos to make columns
    for col, combo in zip(combo_columns, combos):
        table[col] = table[[x for x in combo]].product(axis=1)
    return table


            
"""      
exp.effects_table.loc[:,'a_contrast']   
np.zeros(5)

exp = Experiment2N(['a','b','c',])
exp_table = exp.generate_table()    
exp_table['y'] = [60, 72, 54, 68, 52, 83, 45, 80]
exp_table = exp_table[['a','b','c', 'y']]
effects = exp.calculate_effects(exp_table, 'y')
[''.join(x) for x in exp.contrast_columns]
    

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
'''
test_data = generate_standard_table(3, ['a','b','c'])
for col in test_data.columns:
    test_data[col] = test_data[col].map({0:-1, 1:1})

factors = ['a','b','c']
all_combos = []
for i in range(1, len(factors)+1):
    all_combos.extend(list(combinations(factors,i)))

test_data[[x for x in all_combos[-1]]].prod(axis=1)

test_data['y'] = exp_table.y.values
test_data2 = test_data.copy().drop('y', axis=1)

test_data = test_data.set_index(['a','b','c'])
test_data2 = test_data2.set_index(['a','b','c'])

test_data.merge(test_data2, how='left', left_index=True, right_index=True).reset_index()

test_data  = test_data.reset_index()
test_data.loc[:,'a']*test_data.loc[:,'b']
'''