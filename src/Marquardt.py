#Marquardt.py

from __future__ import print_function, division
import numpy as np
from math import sqrt
from functions import *

EPSILON = 0.00001
MAX_ITERATIONS_NR = 100


def Marquardt(v0, vmin, vmax, x, yObs):
    n = len(v0)
    Lambda0 = 0.01
    vFactor = 2.
    l = np.array([Lambda0]*n)    # damping parameters
    v = np.zeros(n, float)
    for i in range(n): v[i] = v0[i]
    
    nrIter = 1
    maxDiff = 1.0
    sse = norm(v0, x, yObs)
    
    while (maxDiff > EPSILON) and (nrIter < MAX_ITERATIONS_NR):  
        diff = LeastSquares(l, v, vmin, vmax, x, yObs)
        maxDiff = max(abs(diff))
        v_new = computeNewParameters(v, vmin, vmax, diff, l, vFactor)
        sse_new = norm(v_new, x, yObs)
        
        if (sse_new < sse):
            sse = sse_new
            for i in range(n): v[i] = v_new[i]
            l /= vFactor
        else:
            l *= vFactor
        nrIter += 1
    
    yEst = estimate(v, x)
    print ("iterations nr:", nrIter)
    print ("sum of squared residuals:", "{0:.1f}".format(sse))
    print ("a:", "{0:.3f}".format(v[0]))
    print ("b:", "{0:.3f}".format(v[1]))
    print ("c:", "{0:.3f}".format(v[2]))
    return v, yEst


def estimate(v, x):
    bulkDensity = np.zeros(len(x))
    for i in range(len(bulkDensity)):
        bulkDensity[i] = getBulkDensityCurioni(x[i], v[0], v[1], v[2])
        
    return bulkDensity
   

def computeNewParameters(v, vmin, vmax, diff, l, factor):
    n = len(v)
    v_new = np.zeros(n, float)
    for i in range(n):
        v_new[i] = v[i] + diff[i]
        if (v_new[i] > vmax[i]):
            v_new[i] = vmax[i]
            l[i] *= factor
        if (v_new[i] < vmin[i]):
            v_new[i] = vmin[i]
            l[i] *= factor  
    return v_new


def norm(v, x, yObs):
    yEst = estimate(v, x)
    norm = 0
    for i in range(len(x)):
        dy = yObs[i] - yEst[i]
        norm += (dy*dy)
    return norm


def LeastSquares(l, v, vmin, vmax, x, yObs):
    n = len(v)
    m = len(x)
    p = np.resize(np.zeros(n, float),(n,m))
    a = np.resize(np.zeros(n, float),(n,n))
    z = np.zeros(n, float)
    g = np.zeros(n, float)
    v1 = np.zeros(n, float)
    diff = np.zeros(n, float)
    
    for i in range(n): v1[i] = v[i]
    est = estimate(v, x)                           
       
    for i in range(n):
        change = (vmax[i] - vmin[i]) * 0.01
        v1[i] +=  change  
        # get a new set of estimates                              
        yEst = estimate(v1, x)                      
        v1[i] -= change                                 
        for j in range(m):
            # compute derivatives
            if (change > 0):
                p[i][j] = (yEst[j] - est[j]) / change
            else:
                p[i][j] = 0       
    
    for i in range(n):
        for j in range(i, n):
            a[i][j] = 0
            for k in range(m):
                a[i][j] = a[i][j] + p[i][k] * p[j][k]
        z[i] = sqrt(a[i][i]) + EPSILON

    for i in range(n):
        g[i] = 0
        for k in range(m):
            g[i] = g[i] + p[i][k] * (yObs[k] - est[k])
        g[i] = g[i] / z[i]
        for j in range(i, n):
            a[i][j] = a[i][j]  / (z[i] * z[j])
       
    for i in range(n):
        a[i][i] = a[i][i] + l[i]
        for j in range(i+1, n):
            a[j][i] = a[i][j]

    for j in range(n-1):
        pivot = a[j][j]
        for i in range(j+1, n):
            mult = a[i][j] / pivot
            for k in range(j+1, n): a[i][k] -= mult * a[j][k]
            g[i] -=  mult * g[j]

    diff[n-1] = g[n-1] / a[n-1][n-1]
     
    for i in range(n-2, -1, -1):
        top = g[i]
        for k in range(i+1, n):
            top -= a[i][k] * diff[k]
        diff[i] = top / a[i][i]
    
    for i in range(n): 
        diff[i] /=  z[i]

    return diff

