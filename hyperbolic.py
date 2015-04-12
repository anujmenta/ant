
import numpy as np

def funct(c):
    dt = 0.04
    dx = 0.05
    nu = c*dt/dx

    # Max values of X and Y are 4
    n = int(40/dt)
    m = int(4/dx)

    X = np.linspace(-2, 2, m)
    T = np.linspace(-5, 5, n)

    U = np.zeros((m, n), dtype=np.float32)  # U at nth iteration
    # U at (n + 1)th iteration, yet to be calculated


    # Intializing U
    # u = 0; for x = 0; y = 0
    # u = x**2 along y = 4; 0< x < 4
    # u = 16*y along x = 4; 0 < y < 4
    # U[0:int(1/dx),0] = 0
    U[0:int(2/dx), 0] = X[0:int(2/dx)]+2
    U[int(2/dx):, 0] = 2 - X[int(2/dx):]
    # U[int(1/dx):,0] = 4-X[int(1/dx):]**2
    # U[int(5/dx):,0] = 0

    for i in range(0, n-1):
        for j in range(0, m):
            U[j, i+1] = U[j,i] - nu*(U[j,i] - U[(j-1+m)%m, i])

    return U



def lax_wend(c):
    dt = 0.04
    dx = 0.05
    nu = c*dt/dx

    # Max values of X and Y are 4
    n = int(40/dt)
    m = int(4/dx)

    X = np.linspace(-2, 2, m)
    T = np.linspace(-5, 5, n)

    U = np.zeros((m, n), dtype=np.float32)  # U at nth iteration
    # U at (n + 1)th iteration, yet to be calculated


    # Intializing U
    # u = 0; for x = 0; y = 0
    # u = x**2 along y = 4; 0< x < 4
    # u = 16*y along x = 4; 0 < y < 4
    # U[0:int(1/dx),0] = 0
    U[0:int(2/dx), 0] = X[0:int(2/dx)]+2
    U[int(2/dx):, 0] = 2 - X[int(2/dx):]
    # U[int(1/dx):,0] = 4-X[int(1/dx):]**2
    # U[int(5/dx):,0] = 0

    for i in range(0, n-1):
        for j in range(0, m):
            U[j, i+1] = U[j,i] - nu*(U[(j+1)%m,i] - U[(j-1+m)%m, i])/2 + (nu**2)*(U[(j+1)%m, i] - 2*U[j,i] + U[(j-1+m)%m, i])/2

    return U
