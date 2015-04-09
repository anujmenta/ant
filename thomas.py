# Author: Harsh Gupta
# Roll No: 12MA20017

# Using Thomas algorithm to solve the linear BVP
import sys
import numpy as np
import matplotlib.pyplot as plt


def read_BVP(infile):
    """
    Read the input file to returns a dict
    The equation is in form: A y" + B y' + C y - D = 0
    and the boundary_conditions are y(l) = y_l and y(r) = y_r

    """
    coeffs = [int(x) for x in infile.readline()[:-1].split(' ')]
    initial_cond = [int(x) for x in infile.readline()[:-1].split(' ')]
    final_cond = [int(x) for x in infile.readline()[:-1].split(' ')]
    print("Enter h")
    h = float(infile.readline()[:-1])
    return (coeffs, initial_cond, final_cond, h)


def thomas_algorithm(a, b, c, d):
    """
    Solves the Tridiagonal Linear System
          --             -- -- --   -- --
          |b_1 c_1        | |f_1|   |d_1|
          |a_2 b_2 c_2    | | . |   | . |
          |    a_3 . . .  | | . | = | . |
          |               | |   |   |   |
          |               | |   |   |   |
          |       a_n b_n | |f_n|   |d_n|
          --             -- -- --   -- --
    """
    assert len(a) == len(b) == len(c) == len(d)
    N = len(c)
    c_ = [None for i in range(N)]
    d_ = [None for i in range(N)]
    f = [None for i in range(N)]
    c_[0] = c[0]/b[0]
    d_[0] = d[0]/b[0]

    for i in range(1, N):
        c_[i] = c[i]/(b[i] - a[i]*c_[i-1])
        d_[i] = (d[i] - a[i]*d_[i-1])/(b[i] - a[i]*c_[i-1])

    f[N-1] = d_[N-1]
    for i in range(N-2, -1, -1):
        f[i] = d_[i] - c_[i]*f[i+1]

    return f


def solve_BVP(coeffs, initial_cond, final_cond, h):
    A, B, C, D = coeffs
    l, y_l = initial_cond
    r, y_r = final_cond

    N = (r - l)/h
    assert int(N) == N
    N = int(N)

    a = [None for i in range(N-1)]
    b = [None for i in range(N-1)]
    c = [None for i in range(N-1)]
    d = [None for i in range(N-1)]

    for i in range(N-1):
        a[i] = A/(h**2) - B/(2*h)
        b[i] = -2*A/(h**2) + C
        c[i] = A/(h**2) + B/(2*h)
        d[i] = D
        if i == 0:
            d[i] = D - y_l*(A/(h**2) - B/(2*h))

    d[N-2] = D - y_r*(A/(h**2) + B/(2*h))

    f = [y_l] + thomas_algorithm(a, b, c, d) + [y_r]
    return f


if __name__ == "__main__":
    input_file = open("in", "r")
    coeffs, initial_cond, final_cond, h = read_BVP(input_file)
    f = solve_BVP(coeffs, initial_cond, final_cond, h)
    print(f)
