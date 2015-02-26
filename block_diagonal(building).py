import numpy as np

def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    '''
    
    n = len(b)
    #Changed to floats

    aa, bb, cc, dd = [0]*n, [0]*n, [0]*n, [0]*n

    #initialized

    cc[0] = (np.linalg.inv(b[0]))*c[0]
    dd[0] = (np.linalg.inv(b[0]))*d[0]

    for j in range(1,n):
        bb[j] = b[j] - a[j]*cc[j-1]
        cc[j] = (np.linalg.inv(bb[j]))*c[j]
        dd[j] = (np.linalg.inv(bb[j]))*(d[j]-a[j]*dd[j-1])

    print bb, cc, dd

    for num in range(n-2,-1,-1):
        dd[num] = dd[num] - cc[num]*dd[num+1]
    

    return dd
