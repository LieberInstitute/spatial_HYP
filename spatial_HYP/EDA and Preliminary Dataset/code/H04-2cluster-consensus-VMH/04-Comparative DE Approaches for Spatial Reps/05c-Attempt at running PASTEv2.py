#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  1 16:42:17 2023

@author: bmulvey
"""
import math
import time
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import style
#import paste as pst
import ot
import torch

### manually define functions from paste 2 scripts as available on github master branch 07/01/2023 ...

slice1=sc.read("/Users/bmulvey/Desktop/KM Lab/local_hthspatial/analysis/data/H04-2cluster-consensus-VMH/Comparative DE approaches/d1s1.h5ad")
slice2=sc.read("/Users/bmulvey/Desktop/KM Lab/local_hthspatial/analysis/data/H04-2cluster-consensus-VMH/Comparative DE approaches/d1s2.h5ad")

slice1.obsm['spatial']=np.genfromtxt("/Users/bmulvey/Desktop/KM Lab/local_hthspatial/analysis/data/H04-2cluster-consensus-VMH/Comparative DE approaches/d1s1_coord.csv",delimiter=',')

slice2.obsm['spatial']=np.genfromtxt("/Users/bmulvey/Desktop/KM Lab/local_hthspatial/analysis/data/H04-2cluster-consensus-VMH/Comparative DE approaches/d1s2_coord.csv",delimiter=',')


pi12 = partial_pairwise_align(slice1,slice2,s=0.9)

### hangs at:
# Traceback (most recent call last):

#   Cell In[46], line 1
#     pi12 = partial_pairwise_align(slice1,slice2,s=0.9)

#   Cell In[42], line 283 in partial_pairwise_align
#     pi, log = partial_fused_gromov_wasserstein(M, D_A, D_B, a, b, alpha=alpha, m=m, G0=G_init, loss_fun='square_loss', armijo=armijo, log=True, verbose=verbose)

#   Cell In[42], line 158 in partial_fused_gromov_wasserstein
#     gamma = ot.optim.solve_1d_linesearch_quad(a, b, c)

# TypeError: solve_1d_linesearch_quad() takes 2 positional arguments but 3 were given


# armijo
# Traceback (most recent call last):

#   Cell In[47], line 1
#     armijo

# NameError: name 'armijo' is not defined    



