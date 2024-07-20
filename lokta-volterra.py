from copy import deepcopy
import numpy as np
import time
import matplotlib.pyplot as plt
import math 



def EulerLoktaVolterra2(courant,pred, a, b, c, d,h):
    """implémentation equation euler implicite 

    Args:
        courant (matrix(1,2)): xn_1
        pred (matrix(1,2)): xn
        a (float): alpha
        b (float): beta
        c (float): delta
        d (float): gamma
        h (float): pas de temps

    Returns:
        matrix((1,2)): équations au pas de temps xn_1 (en fonction de xn et des autres paramètres)
        
    """
    
    xn = pred[0]
    yn = pred[1]
    
    xn_1 = courant[0]
    yn_1 = courant[1] 
    
    return np.array([xn_1-xn-h*a*xn_1+h*b*xn_1*yn_1,
                     yn_1-yn-h*c*xn_1*yn_1+h*d*yn_1])
    
def JacobienLoktaVolterra2(courant,a,b,c,d,h):
    
    return np.array([[1-h*a+h*b*courant[1],h*b*courant[0]],
                     [-h*c*courant[1],1-h*c*courant[0]+h*d]])

def NewtonNDim(EulerLoktaVolterra2,JacobienLoktaVolterra2,z0,eps,Nmax,a,b,c,d,h):
    
    x = z0
    
    fx = EulerLoktaVolterra2(x,z0,a,b,c,d,h)
    
    J = JacobienLoktaVolterra2(x,a,b,c,d,h)
    
    k = 0

    while(np.linalg.norm(fx) > eps and k<Nmax ):
        
        k = k + 1 
        
        x = x - np.linalg.solve(J, fx)
        
        fx = EulerLoktaVolterra2(x,z0,a,b,c,d,h)
        
        J = JacobienLoktaVolterra2(x,a,b,c,d,h)
        
    return x


a = 5
b = 0.1
c = 0.01
d = 5
h = 0.1
Nmax = 1000
u = np.array([10,4])



temp=[]
xSol=[]
ySol=[]
sol=u

for t in range(0,1000):
    
    temp.append(t)
    sol = NewtonNDim(EulerLoktaVolterra2,JacobienLoktaVolterra2,sol,1e-7,Nmax,a,b,c,d,h)
    xSol.append(sol[0])
    ySol.append(sol[1])

plt.figure("Lokta-Volterra 2 équations",figsize=(8.0,6.4))
plt.plot(temp,xSol,label="les proies")
plt.plot(temp,ySol,label="les prédateurs")
plt.axhline(y=a/b, color='r', linestyle='--',label="équilibre proie")
plt.axhline(y=d/c, color='g', linestyle='--',label="équilibre prédateur")
plt.xlabel("Temps")
plt.ylabel("Effectif des populations en nombres d'individus")
plt.grid()
plt.title("Lokta-Volterra 2 équations - Euler Implicite")
plt.legend()
plt.show()




