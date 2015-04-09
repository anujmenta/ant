#Anuj Menta
#12MA20007

#Input the arrays A, B, C, D to the TDMASolver function
#Example: y'' = x+y h = 0.25
# Obtains the array A = [0, 16, 16], B = [-33, -33, -33], C = [16, 16, 0], D = [0.25, 0.5, 0.75] (Inputs)
# TDMAsolver(A,B,C,D) --> Gives the result Y = [-0.035, -0.056, -0.05]

import numpy as np
import matplotlib.pyplot as plt
from math import exp
import time

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

    #Thomas Algorithm
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

    #Thomas end (back substitution)

    return dd

#Input co-efficients of the equation
a_inp = raw_input("Enter the coeffiecient of y' : ")
b_inp = raw_input("Enter the coeffiecient of y : ")
c_inp = raw_input("Enter the coeffiecient of -1 : ")

a1, in_value = raw_input("Enter the intial values of x and y in the format x, y(x): ").split()
b1, b_value = raw_input("Enter the boundary values of x and y in the format x, y(x): ").split()

a1 = float(a1)
b1 = float(b1)

#Functions for evaluating the values of co-efficients
def a(x):
    return float(eval(str(eval(a_inp))))
def b(x):
    return float(eval(str(eval(b_inp))))
def c(x):
    return float(eval(str(eval(c_inp))))

#The co-effiecients of the equation in the discretized form
def A(h, x):
    return round(((1/h**2)+(a(x))/(2*h)),2)
def B(h):
    return round(((b(h))-(2/h**2)),2)
def C(h, x):
    return round(((1/h**2)-(a(x))/(2*h)),2)

#Number of iterations?
n1 = int(raw_input("Start Value of n : "))
n2 = int(raw_input("End Value of n : "))

#For loop for plotting the graph
for n_iter in range(n1,n2+1):

    h=1/float(n_iter) #standard notation (h)

    #Generating the arrays (Tridiagonal Matrix)
    A_arr = [A(h, x) for x in np.arange(a1+h,b1,h)]
    B_arr = [B(h)]*(int(((b1-a1)/h)-1))
    C_arr = [C(h, x) for x in np.arange(a1+h,b1,h)]
    D_arr = [c(round(x,3)) for x in np.arange(a1+h,b1,h)]
    D_arr[0]+= C_arr[0]*(-1.)*(float(in_value))
    #print C_arr[0]*(-1.)*(float(in_value))
    D_arr[-1]+= A_arr[0]*(-1.)*(float(b_value))
    #print A_arr[0]*(-1.)*(float(b_value))
    #D_arr[-1] += 16

    #For some cases like n=7, i.e h=1/7 it so happens that number of
    #elements goes haphazrd
    if (len(A_arr)!=len(D_arr)):
        print n_iter
    #Solve the Tri-Diagonal Matrix
    D_sol = TDMAsolver(A_arr, B_arr, C_arr, D_arr)
    h_sol = [round(i,3) for i in np.arange(a1+h, b1, h)]

    #Plotting the graph
    try:
        plt.plot(h_sol,D_sol)
    except:
        pass
    time.sleep(0.1)
    plt.pause(0.0001)

#Naming the graph and properly labeling it.
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.title("Graphical comparision for the equation y''+("+a_inp+")y'+("+b_inp+")y = "+c_inp+" \nfor n ="+str(n1)+" to n ="+str(n2)+". Here h = (1/n)")
plt.savefig('graph_tdma.pdf', format = 'pdf')

#End of program.
