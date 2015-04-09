# Solving Poisson's equation
# u_xx + u_yy = 0
# For the boundary conditions
# u = 0; for x = 0; y = 0
# u = x**2 along y = 4; 0< x < 4
# u = 16*y along x = 4; 0 < y < 4

import numpy as np

dx = 0.05
dy = 0.05

# Max values of X and Y are 4
n = int(4/dx)
m = int(4/dy)
X = np.linspace(0, 4, n)
Y = np.linspace(0, 4, m)

U = np.zeros((n, m), dtype=np.float32)  # U at nth iteration
# U at (n + 1)th iteration, yet to be calculated


# Intializing U
U[n-1, :] = X**2
U[:, m-1] = 16*Y


U_ = np.zeros((n, m), dtype=np.float32)
for k in range(40):  # number of iterations
    U_ = np.zeros((n, m), dtype=np.float32)
    U_[n-1, :] = X**2
    U_[:, m-1] = 16*Y
    for i in range(1, n - 1):
        for j in range(1, m - 1):
            U_[i, j] = (1/(2*(dx**2 + dy**2)))*((dy**2)*(U_[i-1, j] + U[i+1, j]) + (dx**2)*(U_[i, j-1] + U[i, j+1]))
    U = U_.copy()
