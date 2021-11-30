import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate.odepack import odeint

##### problem i #####

#generating vector field mesh
x,y = np.meshgrid(np.linspace(-10,250,20), np.linspace(-10,130,20))

Lambda = 32.221
b = 1.3847
Mu = 27.432
Mu2 = Mu+5
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

k1,x1,y1 = find_k(f1[0],f1[1],Lambda,b,Mu,c)
k2,x2,y2 = find_k(f2[0],f2[1],Lambda,b,Mu,c)

equation1 = ((U-x1)**2)/(b/(c*Lambda)) + ((V-y1)**2)/(c/(b*Mu))
equation2 = ((U-x2)**2)/(b/(c*Lambda)) + ((V-y2)**2)/(c/(b*Mu))
plt.contour(U,V,equation1, [k1])
plt.contour(U,V,equation2, [k2])
plt.show()

##### problem ii #####

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

##### problem iii #####

plt.contour(U,V,equation1, [k1])
plt.contour(U,V,equation2, [k2])
plt.plot(sol[:,0], sol[:,1], label='[21,44]')
plt.plot(sol1[:,0], sol1[:,1], label='[72,23]')
plt.plot(sol2[:,0], sol2[:,1], label='[2,16]')
plt.xlabel("Hares")
plt.ylabel("Lynx")
plt.legend()
plt.show()

##### problem vi #####

sol3 = odeint(LotkaVolterra, f1, t, args=(Lambda, b, Mu2, c))
sol4 = odeint(LotkaVolterra, f2, t, args=(Lambda, b, Mu2, c))
sol5 = odeint(LotkaVolterra, f3, t, args=(Lambda, b, Mu2, c))

x,y = np.meshgrid(np.linspace(-10,300,20), np.linspace(-10,130,20))
dx = x*(Lambda-b*y)
dy = y*(-Mu+c*x)
t = np.linspace(0,1,1000)

n = np.sqrt(dx**2 + dy**2)
u,v = dx/n, dy/n

plt.plot(sol3[:,0], sol3[:,1],label='[21,44]')
plt.plot(sol4[:,0], sol4[:,1],label='[72,23]')
plt.plot(sol5[:,0], sol5[:,1],label='[2,16]')
plt.xlabel("Hares")
plt.ylabel("Lynx")
plt.legend()
plt.quiver(x,y,u,v)
plt.show()

h1 = Lambda*math.log(f1[1])+ Mu*math.log(f1[0]) - b*f1[1] - c*f1[0]
h2 = Lambda*math.log(f2[1])+ Mu*math.log(f2[0]) - b*f2[1] - c*f2[0]
h3 = Lambda*math.log(f3[1])+ Mu*math.log(f3[0]) - b*f3[1] - c*f3[0]



x1,y1 = np.meshgrid(np.linspace(0,260,1000), np.linspace(0,110,1000))
f = Lambda*np.log(y1) + Mu*np.log(x1) - b*y1 - c*x1

plt.contour(x1,y1,f,[h1],label='[21,44]',colors='blue')
plt.contour(x1,y1,f,[h2],label='[72,23]',colors='orange')
plt.contour(x1,y1,f,[h3],label='[2,16]',colors='green')
plt.xlabel("Hares")
plt.ylabel("Lynx")
plt.legend()
plt.show()

##### problem vii #####

u_0,v_0 = -10,10
K = 400
q = 0.88

#equilibrium points
x_star = (Mu/c)**(1/q)
y_star = (Lambda/b)*((Mu/c)**((1/q)-1))*(1-(1/K)*x_star)
#y_star = (Lambda*((Mu/c)**(((1-q)/q)))-((1/K)*((Mu/c)**(((2-q/q))))))/b 


def solve_linear_sys(s,t):
    u,v = s[0],s[1]
    uu = Lambda*(1-(2*x_star)/K)-b*y_star*q*(x_star**(q-1))
    uv = -b*(x_star**q)
    vu = c*y_star*q*(x_star**(q-1))
    vv = -Mu + c*(x_star**q)
    dudt = uu*u + uv*v
    dvdt = vu*u + vv*v
    return [dudt,dvdt]

t = np.linspace(0,3,1000)
new_sol = odeint(solve_linear_sys, [u_0,v_0], t)

plt.plot(new_sol[:,0] + x_star, new_sol[:,1] + y_star)
plt.xlabel("Hares")
plt.ylabel("Lynx")
plt.legend()
plt.show()



##### problem viii #####






##### problem ix #####






##### problem x #####






##### problem xi #####





