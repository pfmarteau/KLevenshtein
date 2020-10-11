#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 16:46:45 2020
@author: pfm
"""
from scipy.spatial.distance import cdist
import numpy as np
# Filename: kLevenshtein.py
# Python source code for the "Kernelized" Levenshtein distance similarity (Kernelization as defined in the reference below).
# Author: Pierre-Francois Marteau
# Version: V1.0 du 27/09/2020, 
# Licence: GPL
# ******************************************************************
# This software and description is free delivered "AS IS" with no 
# guaranties for work at all. Its up to you testing it modify it as 
# you like, but no help could be expected from me due to lag of time 
# at the moment. I will answer short relevant questions and help as 
# my time allow it. I have tested it played with it and found no 
# problems in stability or malfunctions so far. 
# Have fun.
# *****************************************************************
# Please cite as:
# @article{marteau:hal-00486916,
#   AUTHOR = {Marteau, Pierre-Francois and Gibet, Sylvie},
#   TITLE = {{On Recursive Edit Distance Kernels with Application to Time Series Classification}},
#   JOURNAL = {{IEEE Transactions on Neural Networks and Learning Systems}},
#   PAGES = {1-14},
#   YEAR = {2014},
#   MONTH = Jun,
#   KEYWORDS = {Elastic distance, Time warp kernel, Time warp inner product, Edit distance kernel, Definiteness, Time series classification, SVM},
#   DOI = {10.1109/TNNLS.2014.2333876},
#   URL = {http://hal.inria.fr/hal-00486916}
# } 
# 
'''
# input A: first ordinal time series: array of array (nx1), n is the number of characters translated as ordinal
# intput B: second time series: array of array (nx1), n is the number of characters translated as ordinal
# input local_kernel: matrix of local kernel evaluations
'''
def kLevenshtein_lk(A, B, local_kernel):
    d=np.shape(A)[1]
    Z=[np.zeros(d)]
    A = np.concatenate((Z,A), axis=0)
    B = np.concatenate((Z,B), axis=0)
    [la,d] = np.shape(A)
    [lb,d] = np.shape(B)
    DP = np.zeros((la,lb))
    DP1 = np.zeros((la,lb));
    DP2 = np.zeros(max(la,lb));
    l=min(la,lb);
    DP2[1]=1.0;
    for i in range(1,l):
        DP2[i] = local_kernel[i-1,i-1];

    DP[0,0] = 1;
    DP1[0,0] = 1;
    n = len(A);
    m = len(B);

    for i in range(1,n):
        DP[i,1] = DP[i-1,1]*local_kernel[i-1,2];
        DP1[i,1] = DP1[i-1,1]*DP2[i];

    for j in range(1,m):
        DP[1,j] = DP[1,j-1]*local_kernel[2,j-1];
        DP1[1,j] = DP1[1,j-1]*DP2[j];

    for i in range(1,n):
        for j in range(1,m): 
            lcost=local_kernel[i-1,j-1];
            DP[i,j] = (DP[i-1,j] + DP[i,j-1] + DP[i-1,j-1]) * lcost;
            if i == j:
                DP1[i,j] = DP1[i-1,j-1] * lcost + DP1[i-1,j] * DP2[i] + DP1[i,j-1]  *DP2[j]
            else:
                DP1[i,j] = DP1[i-1,j] * DP2[i] + DP1[i,j-1] * DP2[j];
    DP = DP + DP1;
    return DP[n-1,m-1]




''''
# function kLevenshteind(A, B, sigma, epsilon)
# Dynamic programming implementation of kLevenshtein kernel
# input A: first multivariate time series: array of array (nxd), n is the number of sample, d is the dimension of each sample
# intput B: second multivariate time series: array of array (nxd), n is the number of sample, d is the dimension of each sample
# input sigma: >0, used in the exponential local kernel 
# input epsilon: 1 > epsilon > 0, used in the exponential local kernel 
# output similarity: similarity between A and B (the higher, the more similar)
'''
def kLevenshtein(A, B, sigma = .1, epsilon = 1e-3):
    A=str2array(A)
    B=str2array(B)
    distance = cdist(A, B, 'hamming')
    local_kernel = (np.exp(-distance/sigma)+epsilon)/(3*(1+epsilon))
    return kLevenshtein_lk(A,B,local_kernel)

def str2array(s):
    out=[]
    for e in s:
        out.append([ord(e)])
    return np.array(out)

# Simple test
''' This example shows that it is possible by selecting a high epsilon meta-parameter to force the string kernel to evaluate not only 
the global alignments of the two compared strings but also the alignments of all matching substrings. Hence, for Klevenshtein the maximum similarity 
is reached for  the strings *[the little big man]* and *[the big little man]*, while for the Levenshtein distance, the minimal distance is reached for the strings
*[the little big man]* and *["the little mairmaid"]*. This demonstrates that the KLevenshtein kernel can be more resilinet than the Levenshtein's distance to 
string permutations.
'''
if __name__ == '__main__':
    import Levenshtein
    
    A = "the little mairmaid"
    B = "the little big man"
    C = "the big little man"
    D = "an old mairmaid at sea"
    sigma=5
    epsilon=1
    L=[A,B,C,D]
    GLev=np.zeros((len(L),len(L)))
    GKLev=np.zeros((len(L),len(L)))
    for i in range(len(L)):
        for j in range(i+1, len(L)):
            GLev[i,j]=Levenshtein.distance(L[i],L[j])
            GKLev[i,j]=kLevenshtein(L[i],L[j],sigma=sigma,epsilon=epsilon)
            print("[%s] v.s. [%s] | lev= %d | klev= %1.4f" % (L[i], L[j], GLev[i,j], GKLev[i,j]))
