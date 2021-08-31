#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 13:36:05 2020

@author: seifhejazine
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
matplotlib.rcParams['figure.dpi'] = 300


#list of potentials to trial with, commment out what you are not using
#if you pick the final potential make sure to set dE<0
#def V(x):
#    return 0.5*(x**2) - 0.1*x**3
#def V(x): #morse potential with all constants = 1
 #   return 2*(np.exp(-2*x)-2*np.exp(-1*x))
l=5 #Pöschl–Teller potential
def V(x):
    return l*(l+1)/(-2*(np.cosh(x))**2)

sg = True #True to plot on SameGraph
h = 0.001 #x differential
xout = np.arange(-6.0, 6.0, h)

def f(E, x): #k(x) in numerov's method
    return 2* (E - V(x))

def normal(x): #normalisation function
    norm = np.linalg.norm(x)
    return x/norm    

def numerov(E):
    yout = np.zeros(12000)
    N  = 6000
    y = 0.0
    x  =  -1*(N-2)*h
    b = (h**2)/12 #just a factor
                                                                         
    
    kd2 = f(E, x-2*h)
    kd1 = f(E, x-h)
    yd2 = 0.0
    first_guess = 0.01
    yd1 = first_guess
    
    n = -1*N+2

###Numerov
    m = 0
    while n<N+2:
        n += 1
        x += h
        kx = f(E, x)
        y = (2.0*yd1*(1.0-(5.0*b*kd1)) - yd2*(1.0+b*kd2))/(1+b*kx)
        
        yout[m] = y
        m += 1
        
        yd2 = yd1
        yd1 = y
        kd2 = kd1
        kd1 = kx
        
    return normal(yout)


def bisect(El, Eu, n):
    a = El
    b = Eu
    itr = 1
    while itr<= maxm:
        c = (a+b)/2.0
        psi0 = numerov(c)
        if abs(psi0[-1])**2 < tol:
            print('possible match')
            return c
            break
        itr += 1
        print('itr= ', itr)
        if itr == maxm:
            print('reached max iterations')
            return c
            break
        psi1 = numerov(a)
        if np.sign(psi0[-1]) == np.sign(psi1[-1]):
            a = c
        else:
            b = c

tol = 1.0e-8 #tolerance for value to be considered 0
maxm = 200 #max iterations of bisection

dE = 0.1 #energy differential
E0 = 0.01 #initial energy

M = 6
cstr=["#332288", "#88CCEE", "#44AA99", '#117733', '#999933', "#DDCC77", "#CC6677", '#882255', "#AA4499"]


def main():
    elist=[]
    wvlist=[]
    E1 = E0
    n = 1
    while n <= M:
        E2 = E1 + dE
        print('E1: ', E1,  ' E2: ',  E2)
        psi1 = numerov(E1)
        psi2 = numerov(E2)
        if np.sign(psi1[-1]) == np.sign(psi2[-1]):
            E1 = E2
            print('none in this interval')
        else:
            print('beginning bisection')
            E = bisect(E1, E2, n)
            psi1 = numerov(E)
            wvlist.append(psi1)
            elist.append(E)
            print(E)
            E1 = E2
            n += 1
    return wvlist, elist #returns list of wavefunctions and list of energyvalues

def plotting(list): #plots all psi on same axis
    n=0
    figure(num=None, figsize=(8, 6), dpi=380, facecolor='w', edgecolor='k')
    plt.title('QHO with 3rd order perturbation, max_n=6')  
    for i in list[0]:
        plt.plot(xout,i, color=cstr[n%9])
        n+=1
    print(list[1])
    plt.show()

            
def plot(E): #plots the wavefunction for a given energy 
      psi1=numerov(E)
      plt.plot(xout,psi1)
      plt.show()
      return psi1
  
def plotall(list): #plots each psi individually
    n=0
    for i in list[0]:
        figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
        plt.plot(xout,i, color=cstr[n%9])
        plt.show
        n+=1


if sg:
    plotting(main())
else:
    plotall(main())
    
