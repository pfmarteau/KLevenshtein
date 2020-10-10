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
#   KEYWORDS = {Elastic distance, Time warp kernel, Time warp inner product, Definiteness, Time series classification, SVM},
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
if __name__ == '__main__':
    import Levenshtein
    
    A = "the little mairmaid"
    B = "the old man and the sea"
    C = "an old mairmaid at sea"

    print("kLevenshtein(A,B)=", kLevenshtein(A,B,.1))
    print("kLevenshtein(A,C)=", kLevenshtein(A,C,.1))
    print("kLevenshtein(B,C)=", kLevenshtein(B,C,.1))
    
    print("Levenshtein(A,B)=", Levenshtein.distance(A,B))
    print("Levenshtein(A,C)=", Levenshtein.distance(A,C))
    print("Levenshtein(B,C)=", Levenshtein.distance(B,C))
