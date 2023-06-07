# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 21:55:59 2023

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

mu = 4*np.pi*10**(-7)
s = 10
L = 10
c = 5
r = 0.05
I_0 = 1
t = 0
dt = 0.03
Nt = 100
N = 100

def E(t):
    return -10**7*mu*I_0*c/(2*np.pi*((c*t)**2 - s**2)**0.5)

pos_0 = np.zeros(shape=(2*N+1, 2))
pos = np.zeros(shape=(2*N+1, 2))
vec = np.zeros(shape=(2*N+1, 2))
circle = [] #circles

for i in range(2*N+1):
    if i <= N:
        pos_0[i][1] = i/L
        pos[i][1] = i/L
    else:
        pos_0[i][1] = -(i-N)/L
        pos[i][1] = -(i-N)/L
    circle.append(plt.Circle((pos[i][0], pos[i][1]), r))
    vec[i][0] = s/(pos[i][1]**2 + s**2)**0.5
    vec[i][1] = -pos[i][1]/(pos[i][1]**2 + s**2)**0.5

fig, axes = plt.subplots()
axes.set_aspect(1)
axes.set_xlim([-2.5, 13.5])
axes.set_ylim([-L, L])
axes.plot([0, 0], [-L, L], color='black')

plt.pause(6)

for i in range(2*N+1):
    axes.add_artist(circle[i])

for i in range(Nt):
    t = t + dt
    plt.cla()
    axes.set_aspect(1)
    axes.set_xlim([-2.5, 13.5])
    axes.set_ylim([-L, L])
    axes.plot([0, 0], [-L, L], color='black')
    for j in range(2*N+1):
        if pos[j][0] < s:
            pos[j][0] = pos[j][0] + c*vec[j][0]*dt
            pos[j][1] = pos[j][1] + c*vec[j][1]*dt
        circle[j] = plt.Circle((pos[j][0], pos[j][1]), r, color='red')
        axes.add_artist(circle[j])
        axes.plot([pos_0[j][0], pos[j][0]], [pos_0[j][1], pos[j][1]], color='blue', linewidth=0.1, zorder=0)
    if t > 2:
        axes.arrow(10, 0, 0, E(t).real, color='black', width=0.05)
    plt.show()
    plt.pause(0.01)
        
        