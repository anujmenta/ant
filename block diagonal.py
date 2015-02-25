import numpy as np

def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    '''
    
    n = len(b)
    x = ['a', 'b','c', 'd']
    a = (np.array(a)).astype(float)
    b = (np.array(b)).astype(float)
    c = (np.array(c)).astype(float)
    d = (np.array(d)).astype(float)
    #Changed to floats

    aa, bb, cc, dd = [0]*n, [0]*n, [0]*n, [0]*n

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
