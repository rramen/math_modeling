import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate.odepack import odeint

##### problem 1 #####

x,y = np.meshgrid(np.linspace(-10,250,20), np.linspace(-10,130,20))

Lambda = 32.221
b = 1.3847
Mu = 27.432
c = 0.5433

f1 = [21,44]
f2 = [72,23]
f3 = [2,16]

def find_k(x,y,Lambda,b,Mu,c):
    dotx = Mu/c
    doty = Lambda/b
    u = x-dotx
    v = y-doty
    k =  (u**2)/(b/(c*Lambda)) + (v**2)/(c/(b*Mu))
    return k, dotx, doty

U,V = np.meshgrid(np.linspace(-10,120,1000),np.linspace(-10,120,1000))

k1,x1,y1 = find_k(f1[0], f1[1], Lambda,b,Mu,c)
k2,x2,y2 = find_k(f2[0], f2[1], Lambda,b,Mu,c)

equation1 = ((U-x1)**2)/(b/(c*Lambda)) + ((V-y1)**2)/(c/(b*Mu))
equation2 = ((U-x2)**2)/(b/(c*Lambda)) + ((V-y2)**2)/(c/(b*Mu))
plt.contour(U,V,equation1, [k1])
plt.contour(U,V,equation2, [k2])
plt.show()

##### problem 2 #####

dx = x*(Lambda-b*y)
dy = y*(-Mu+c*x)
t = np.linspace(0,1,1000)

n = np.sqrt(dx**2 + dy**2)
u,v = dx/n, dy/n

def LotkaVolterra(z,t,a,b,c,d):
    x,y = z
    dydt = [a*x-b*x*y,-c*y+d*x*y]
    return dydt

sol = odeint(LotkaVolterra, f1, t, args=(Lambda, b, Mu, c))
sol1 = odeint(LotkaVolterra, f2, t, args=(Lambda, b, Mu, c))
sol2 = odeint(LotkaVolterra, f3, t, args=(Lambda, b, Mu, c))

plt.plot(sol[:,0], sol[:,1])
plt.plot(sol1[:,0], sol1[:,1])
plt.plot(sol2[:,0], sol2[:,1])
plt.quiver(x,y,u,v)
plt.show()

##### problem 3 #####

plt.contour(U,V,equation1, [k1])
plt.contour(U,V,equation2, [k2])
plt.plot(sol[:,0], sol[:,1], label='[21,44]')
plt.plot(sol1[:,0], sol1[:,1], label='[72,23]')
plt.plot(sol2[:,0], sol2[:,1], label='[2,16]')
#plt.quiver(x,y,u,v)
plt.xlabel("Hares")
plt.ylabel("Lynx")
plt.legend()
plt.show()



