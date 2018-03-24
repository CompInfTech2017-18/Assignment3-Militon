# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 19:10:33 2018

@author: kazy
"""

import numpy as np
import matplotlib.pyplot as plt


class themperture:
    
    def __init__(self, a, T0, T0_border, X, Tem, dx, dt, n, method):
        self.a = a
        self.T0_border = T0_border 
        self.X = X
        self.Tem = Tem
        self.dx = dx
        self.dt = dt
        self.Nx = int(X/dx + 1)
        self.Nt = int(Tem/dt + 1)
        self.method = method
        self.eta = a**2*dt/dx**2
        print (self.eta)
        self.analytical_number_of_members = 100
        self.T0 = T0 
        self.T = np.zeros(shape=(self.Nt, self.Nx))
        self.T[0] = T0_border
        self.T[:,0] = T0_border[0] 
        self.T[:,-1] = T0_border[-1] 
        self.start = 1
        self.n = n

    
    def analytic(self):
        greedX, greedTem = np.meshgrid(np.linspace(0, self.X, self.Nx), np.linspace(0, self.Tem, self.Nt))
        for k in range(int(self.n)):
            self.T += 4.*self.T0/np.pi*np.exp(-(np.pi*self.a/self.X*(2*k+1))**2*greedTem)*np.sin(np.pi*(2*k+1)/self.X*greedX)/(2*k+1)
    
    def exscheme(self, n):
        for i in range(self.start, n + 1):
            self.T[i,1:self.Nx-1] = (1.-2.*self.eta)*self.T[i-1,1:self.Nx-1]+self.eta*(self.T[i-1,:self.Nx-2] + self.T[i-1,2:])
        self.start = n + 1
    
    
    def noex(self, n): 
        alpha = np.zeros(self.Nx - 2)
        beta = np.zeros(self.Nx - 2)
        A = -np.ones(self.Nx - 2)
        C = np.ones(self.Nx - 2)*(2./self.eta + 2.)
        B = -np.ones(self.Nx - 2)

        for j in range(n):
            F = self.T[j, 0:self.Nx-2] + (2./self.eta - 2.)*self.T[j, 1:self.Nx-1] + self.T[j, 2:self.Nx]
            F[0] += self.T[j+1, 0]
            F[-1] += self.T[j+1, -1]
            alpha[0] = -B[0]/C[0]
            beta[0] = F[0]/C[0]
            for i in range(self.Nx - 3):              
                alpha[i + 1] = -B[i]/(A[i]*alpha[i] + C[i])
                beta[i + 1] = (F[i] - A[i]*beta[i])/(A[i]*alpha[i] + C[i])
            self.T[j+1][self.Nx-2] = (F[-1] - A[-1]*beta[-1])/(C[-1] + A[-1]*alpha[-1])
            for i in range(self.Nx - 3, 0, -1):
                self.T[j+1][i] = alpha[i]*self.T[j+1][i+1] + beta[i]
                
    def exnewton(self, n):
        self.Te = 30.
        self.h = 5.
        for i in range(self.start, n + 1):
            self.T[i,1:self.Nx-1] = (1.-2.*self.eta)*self.T[i-1,1:self.Nx-1]+self.eta*(self.T[i-1,:self.Nx-2] + self.T[i-1,2:]) - self.h*dt*(self.T[i-1,1:self.Nx-1] - self.Te)
        self.start = n + 1


               
    def methodic(self, n):
        if self.method == 'Analytical':
            self.analytic()
            
        if self.method == 'Explict':
            self.exscheme(n)          
        
        if self.method == 'NonExplict':
            self.noex(n)
        
        if self.method == 'ExplictNewton':
            self.exnewton(n)
 
    def show_plot(self):
        greedX, greedTem = np.meshgrid(np.linspace(0, self.X, self.Nx), np.linspace(0, self.Tem, self.Nt))
        fig, ax0 = plt.subplots(figsize=(10,6))
        ax0.contourf(greedX, greedTem, self.T, levels = np.linspace(0, self.T0, 11))
        plt.show()  
        
    def run(self, snapshots):
        n = 0
        for n1 in snapshots:
            n = n1 - n
            self.methodic(n)
            self.show_plot()    

T0 = 100.  
(a, X, Tem, dx, dt, n) = (0.2, 1., 1., .01, .0001, 100) 
snapshots = [int(Tem/dt)]
T0_border = np.ones(int(X/dx + 1))*T0
T0_border[0] = T0_border[-1]  = 0

obj = themperture(a, T0, T0_border, X, Tem, dx, dt,n, 'Analytical')
obj.run(snapshots)

objective2 = themperture(a, T0, T0_border, X, Tem, dx, dt,n, 'Explict')
objective2.run(snapshots)

T0_border = T0*np.sin(np.pi*np.linspace(0, X, int(X/dx + 1))/X)
T0_border[0] = T0_border[-1]  = 0 


objective3 = themperture(a, T0, T0_border, X, Tem, dx, dt, n, 'Explict')
objective3.run(snapshots)

T0_border = np.ones(int(X/dx + 1))*T0
T0_border[0] = T0_border[-1]  = 0

objective4 = themperture(a, T0, T0_border, X, Tem, dx, dt, n,'NonExplict')
objective4.run(snapshots)

T0_border = np.ones(int(X/dx + 1))*T0
T0_border[:int((X/dx + 1)/2)] = T0/2
T0_border[0] = T0_border[-1]  = 0

objective5 = themperture(a, T0, T0_border, X, Tem, dx, dt,n, 'NonExplict')
objective5.run(snapshots)

T0 = 100.  
(a, X, Tem, dx, dt) = (0.4, 1., 1., .01, .0001) 
T0_border = np.ones(int(X/dx + 1))*T0
T0_border[0] = T0_border[-1]  = 0

objective6 = themperture(a, T0, T0_border, X, Tem, dx, dt,n, 'ExplictNewton')
objective6.run(snapshots)

T0 = 100.  
(a, X, Tem, dx, dt) = (0.4, 1., 1., .01, .0001) 
T0_border = np.ones(int(X/dx + 1))*T0
T0_border[:int((X/dx + 1)/2)] = T0/2
T0_border[0] = T0_border[-1]  = 0

objective7 = themperture(a, T0, T0_border, X, Tem, dx, dt,n, 'ExplictNewton')
objective7.run(snapshots)