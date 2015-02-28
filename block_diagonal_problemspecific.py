import numpy as np

#Input the matrix arrays a, b, c, d to the TDMAsolver function ( Here a is again an array of numpy matrices.)
#Example: y''''+ 81y = 81x^2 h = 0.24
# Solution: [ z1 = [-0.75516624], y1 = [ 0.08216288]],[ z2 = [-1.6098758 ], y2 = [ 0.14072681]]),
# [ z3 = [-1.91138985], y3 = [ 0.12538318]]


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

#A simple function to square a number
def sq(a):
    return a**2

def A(h):
    return np.matrix([[(1/sq(h)), 0], [1, (-2/sq(h))]])

def B(h):
    return np.matrix([[(-2/sq(h)), 81], [1, (4/sq(h))]])

def C(h):
    return np.matrix([[(1/sq(h), 0)], [0, (-2/sq(h))]])

def D(x):
    return np.matrix([[81*sq(x)],[0]])

