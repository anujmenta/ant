#Anuj Menta
#12MA20007

#Input the arrays A, B, C, D to the TDMASolver function
#Example: y'' = x+y h = 0.25
# Obtains the array A = [0, 16, 16], B = [-33, -33, -33], C = [16, 16, 0], D = [0.25, 0.5, 0.75] (Inputs)
# TDMAsolver(A,B,C,D) --> Gives the result Y = [-0.035, -0.056, -0.05]

import numpy as np

def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    '''
    n = len(b)
    for i in range(n):
        a[i] = float(a[i])
        b[i] = float(b[i])
        c[i] = float(c[i])
        d[i] = float(d[i])

    #Changed to floats

    aa, bb, cc, dd = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)

    #initialized

    cc[0] = c[0]/b[0]
    dd[0] = d[0]/b[0]

    for j in range(1,n):
        bb[j] = b[j] - a[j]*cc[j-1]
        cc[j] = c[j]/bb[j]
        dd[j] = (d[j]-a[j]*dd[j-1])/bb[j]

    for num in range(n-2,-1,-1):
        dd[num] = dd[num] - cc[num]*dd[num+1]
    
    for i in range(len(dd)):
        dd[i] = round(dd[i],3)

    return dd
