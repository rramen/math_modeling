import numpy as np
import matplotlib.pyplot as plt

x,y = np.meshgrid(np.linspace(-5,5,20), np.linspace(-5,5,20))

Lambda = 32.221
b = 1.3847
Mu = 27.432
c = 0.5433

u = x*(Lambda-(b*y))
v = y*(-Mu+(c*x))

equation = (u**2)/(b/(c*Lambda)) + (v**2)/(c/(b*Mu))

plt.quiver(x,y,u,v)

x1,y1 = 21,44
x2,y2 = 72,23
x3,y3 = 2,16


plt.show()