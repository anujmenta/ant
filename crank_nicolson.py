# Solving advection equation u_t + c*u_x = v*u_xx
# u(x, 0) = sin(pi*x) for 0 <= x <= 1
# u(0, t) = sin(1, t) = 0
# dx = 0.05 dt = 0.005 v = {0.01, 0.1, 0.5, 1}

import numpy as np



def crank_nikolson(c, v):
    dx = 0.05
    dy = 0.005
    # Range of X and T
    rX = 1
    rT = 1
    n = int(rX/dx)
    m = int(rT/dy)
    X = np.linspace(0, rX, n)
    Y = np.linspace(0, rT, m)

    U = np.zeros((n, m), dtype=np.float32)  # U at nth iteration
    U[0, :] = 0
    U[n-1, :] = 0
    U[:, 0] = np.sin(np.pi*X)
    for k in range(1, n-1):
        A = [None] + [v/dx**2 + c/(2*dx) for i in range(1, m)]
        B = [-v/dx**2 - 2/dt + c/(2*dx)] + [v/dx**2 + c/(2*dt) for i in range(1, m-1)] + [-v/dx**2 - 2/dt - c/(2*dx)]
        C = [v/dx**2 - c/(2*dx) for i in range(1, m)] + [None]
        D = list(U[k, :m-2]*(v/dx**2 + c/(2*dx)) + U[k, 1:m-1]*()
