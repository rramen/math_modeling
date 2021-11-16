import numpy as np
import matplotlib.pyplot as plt

x,y = np.meshgrid(np.linspace(0,100,20), np.linspace(0,40,20))


Lambda = 32.221
b = 1.3847
Mu = 27.432
c = 0.5433


dx = x*(Lambda-b*y)
dy = y*(-Mu+c*x)

n = np.sqrt(dx**2 + dy**2)
u,v = dx/n, dy/n

plt.quiver(x,y,dx,dy)
plt.show()

x1,y1 = 21,44
x2,y2 = 72,23
x3,y3 = 2,16

equation = (u**2)/(b/(c*Lambda)) + (v**2)/(c/(b*Mu))
equation2 = (-Mu*np.log(x)+c*x-Lambda*np.log(y)+b*y)