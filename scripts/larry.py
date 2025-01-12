#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 7 2024

@author: Palash Sashittal
"""

import pandas as pd
import sys
import argparse
import itertools
import math
import numpy as np
import gurobipy as gp

def main(args):

    df_clone_type_mat = pd.read_csv(f'{args.i}', index_col = 0)
    nprogs = args.k
    threshold = args.t
    tree_flag = args.r
    
    # remove the observed potencies below the threshold
    df_clone_type_mat = df_clone_type_mat[df_clone_type_mat['counts'] >= threshold]
    # remove clones with less than two differentiated cell types
    df_clone_type_mat = df_clone_type_mat[df_clone_type_mat.values[:,1:-1].sum(axis=1) > 1]
    
    
    # clone type mat exclusing counts and undifferentiated cells
    curr_clone_type_mat = df_clone_type_mat.values[:,1:-1]
    curr_clone_counts = df_clone_type_mat['counts'].values

    # ILP
    nrows = curr_clone_type_mat.shape[0]
    ntypes = curr_clone_type_mat.shape[1]
    model = gp.Model('LARRY')
    
    ## character matrix variables
    a = model.addVars(nprogs, ntypes, vtype=gp.GRB.BINARY, name='a')
    b = model.addVars(nrows, ntypes, vtype=gp.GRB.CONTINUOUS, name='b')
    ## selection variables
    x = model.addVars(nrows, nprogs, vtype=gp.GRB.BINARY, name='x')
    ## set inclusion variables
    if tree_flag:
        y = model.addVars(nprogs, nprogs, vtype=gp.GRB.CONTINUOUS, name='y')
        z = model.addVars(nprogs, nprogs, vtype=gp.GRB.CONTINUOUS, name='y')
        
    ## one-hot constraint for x
    for row_idx in range(nrows):
        xsum = gp.LinExpr()
        for prog_idx in range(nprogs):
            xsum += x[row_idx, prog_idx]
        model.addConstr(xsum == 1)
        
    ## constraints between b and a
    for row_idx in range(nrows):
        for type_idx in range(ntypes):
            for prog_idx in range(nprogs):
                model.addConstr(b[row_idx, type_idx] >= x[row_idx, prog_idx] + a[prog_idx, type_idx] - 1)
                model.addConstr(a[prog_idx, type_idx] >= x[row_idx, prog_idx] + b[row_idx, type_idx] - 1)    

    ## setting b = 1 based on the clone type mat
    for row_idx in range(nrows):
        for type_idx in range(ntypes):
            if curr_clone_type_mat[row_idx, type_idx] == 1:
                model.addConstr(b[row_idx, type_idx] == 1)    

    ## set inclusion constraints
    if tree_flag:
        for prog1 in range(nprogs):
            for prog2 in range(nprogs):
                if prog1 == prog2:
                    continue
                for type_idx in range(ntypes):
                    model.addConstr(a[prog1, type_idx] <= a[prog2, type_idx] - y[prog1, prog2] + 1)

        for prog1 in range(nprogs):
            for prog2 in range(nprogs):
                if prog1 >= prog2:
                    continue
                for type_idx in range(ntypes):
                    model.addConstr(a[prog1, type_idx] + a[prog2, type_idx] <= 2 - z[prog1, prog2])

        for prog1 in range(nprogs):
            for prog2 in range(nprogs):
                if prog1 >= prog2:
                    continue
                for type_idx in range(ntypes):
                    model.addConstr(y[prog1, prog2] + y[prog2, prog1] + z[prog1, prog2] >= 1)    
                
    ## set obj
    obj_sum = gp.LinExpr()
    for row_idx in range(nrows):
        for type_idx in range(ntypes):
            if curr_clone_type_mat[row_idx, type_idx] == 0:
                obj_sum += curr_clone_counts[row_idx] * b[row_idx, type_idx]
    
    model.setObjective(obj_sum, gp.GRB.MINIMIZE)
    model.setParam(gp.GRB.Param.Threads, 1)
    model.setParam(gp.GRB.Param.Method, 4)    

    model.optimize()
    
    print(model.getObjective().getValue())
    
    # sola = np.reshape(model.getAttr('x', a).values(), (nprogs, ntypes))
    # solx = np.reshape(model.getAttr('x', x).values(), (nrows, nprogs))
    # solb = np.reshape(model.getAttr('x', b).values(), (nrows, ntypes))
    sola = np.reshape(np.array(list(model.getAttr('x', a).values())), (nprogs, ntypes))
    solx = np.reshape(np.array(list(model.getAttr('x', x).values())), (nrows, nprogs))
    solb = np.reshape(np.array(list(model.getAttr('x', b).values())), (nrows, ntypes))
    
    df_sola = pd.DataFrame(sola, columns = df_clone_type_mat.columns[1:-1], dtype=int)
    df_solx = pd.DataFrame(solx, index=df_clone_type_mat.index, columns = range(nprogs), dtype=int)
    df_solb = pd.DataFrame(solb, index=df_clone_type_mat.index, columns = df_clone_type_mat.columns[1:-1], dtype=int)

    prefix = args.o
    df_sola.to_csv(f'{prefix}_prog.csv')
    df_solb.to_csv(f'{prefix}_labeling.csv')
    df_solx.to_csv(f'{prefix}_assignment.csv')

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='csv file with clone type matrix with counts')
    parser.add_argument('-k', type=int, help='number of progenitors', default = 1)
    parser.add_argument('-t', type=int, help='threshold', default = 0)    
    parser.add_argument('-r', help='tree constraint? [No]', default=False, action='store_true')
    parser.add_argument('-o', type=str, help='output prefix', required=True)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    main(args)