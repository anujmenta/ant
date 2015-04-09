# Using Python 3
from __future__ import print_function

print("Program to solve boundary value problem:")
print("Interval for  solving x is (1,2)")
print("y(1) = 2")
print("y(2) = 4")
print("BVP Solution for y'' + 2yy' = 4 + 4x^3")
print("using quasi linearization technique")
print("delta_x's denominator : ")

delta_x = int(input())

B = [2 + 2*i/float(delta_x) for i in range(delta_x + 1)]
X = [None] + [(1 + i/float(delta_x))**3 for i in range(1, delta_x)]

count = None

A = [[None for i in range(delta_x + 1)] for j in range(delta_x + 1)]
D = [None for i in range(delta_x)]

while count is not 0:
    count = 0
    A[0][0] = (-2)*delta_x**2 + (B[2] - B[0])*(delta_x)
    A[0][1] = delta_x**2 + B[1]*delta_x

    D[0] = 4 + 4*X[1] + B[1]*(B[2] - B[0])*delta_x - \
        (B[0]*(delta_x**2 - B[1]*delta_x))

    for i in range(2, delta_x-1):
        A[i-1][i-2] = delta_x**2 - B[i]*(delta_x)
        A[i-1][i-1] = (-2)*delta_x**2 + (B[i+1] - B[i-1])*(delta_x)
        A[i-1][i] = delta_x**2 + B[i]*(delta_x)
        D[i-1] = 4 + 4*X[i] + B[i]*(B[i+1]-B[i-1])*delta_x
