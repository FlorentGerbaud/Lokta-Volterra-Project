import numpy as np
import matplotlib.pyplot as plt

def LoktaVolterra3(courant,pred,a,b,c,d,e,f,g,h,i,j,dt):
    
    """Cette fonction représente les 3 équations du système mise sous 
       la forme d'Euler Implicite. Les inconnues sont respectivement :
       xn_1, yn_1 et zn_1

    Args:
        courant (matrix(1,3)): vecteur contenant les termes courants ici : xn_1, yn_1 et zn_1
        pred (matrix(1,3)): vecteur contenant les termes précédents ici : xn, yn et zn
        a (float): _description_
        b (float): _description_
        c (float): _description_
        d (float): _description_
        e (float): _description_
        f (float): _description_
        g (float): _description_
        h (float): _description_
        i (float): _description_
        j (float): _description_
        dt (float): pas de discrétisation

    Returns:
        matrix((3,3)): matrice contenant les 3 équations à résoudre.
    """
    
    xn = pred[0]
    yn = pred[1]
    zn = pred[2]
    
    
    xn_1 = courant[0]
    yn_1 = courant[1]
    zn_1 = courant[2]

    
    return np.array([xn_1 - xn - dt*xn_1*(a-b*xn_1-c*yn_1),
                     yn_1 - yn - dt*yn_1*(d-e*yn_1-f*xn_1-g*zn_1),
                     zn_1 - zn - dt*zn_1*(h*yn_1-i*zn_1-j)])
    
    
def Jacobien3(courant,a,b,c,d,e,f,g,h,i,j,dt):
    """ Évaluation du Jacobien du systéme en (xn_1, yn_1 et zn_1)

    Args:
        courant (matrix(1,3)): vecteur contenant les termes courants ici : xn_1, yn_1 et zn_1
        a (float): _description_
        b (float): _description_
        c (float): _description_
        d (float): _description_
        e (float): _description_
        f (float): _description_
        g (float): _description_
        h (float): _description_
        i (float): _description_
        j (float): _description_
        dt (float): pas de discrétisation

    Returns:
        matri((3,3)) : Le Jacobien évalué au point courant
    """
    
    xn_1 = courant[0]
    yn_1 = courant[1]
    zn_1 = courant[2]
    
    return np.array([[1-dt*a+2*dt*b*xn_1+dt*c*yn_1,dt*c*xn_1,0],
                     [dt*f*xn_1,1-(dt*d)+(dt*2*e*yn_1)+(dt*f*xn_1)+(dt*g*zn_1),dt*g*yn_1],
                     [0,dt*zn_1*h,(1-dt*yn_1*h)+(2*dt*i*zn_1)+(dt*j)]])

def NewtonNDim(LoktaVolterra3,Jacobien3,z0,eps,Nmax,a,b,c,d,e,f,g,h,i,j,dt):
    """_summary_

    Args:
        LoktaVolterra3 (function): passage de la fonction LoktaVolterra3 en paramètre
        Jacobien3 (function): passage de la fonction Jacobien3 en paramètre
        z0 (matrix(1,3)): vecteur d'initialisation (xn_1,yn_1,zn_1)
        eps (float): précision d'arrêt
        Nmax (int): maximum d'itérations
        a (float): _description_
        b (float): _description_
        c (float): _description_
        d (float): _description_
        e (float): _description_
        f (float): _description_
        g (float): _description_
        h (float): _description_
        i (float): _description_
        j (float): _description_
        dt (float): pas de discrétisation

    Returns:
        x (float) : estimation de la solution de l'équation.
    """
    
    x = z0
    
    fx = LoktaVolterra3(x,z0,a,b,c,d,e,f,g,h,i,j,dt)
    
    J = Jacobien3(x,a,b,c,d,e,f,g,h,i,j,dt)
    
    k = 0

    while(np.linalg.norm(fx) > eps and k<Nmax ):
        
        k = k + 1 
        
        x = x - np.linalg.solve(J, fx)
        
        fx = LoktaVolterra3(x,z0,a,b,c,d,e,f,g,h,i,j,dt)
        
        J = Jacobien3(x,a,b,c,d,e,f,g,h,i,j,dt)
        
    return x



#####################################################
#                                                   #
#                          MAIN                     #
#                                                   #
#####################################################

# _________________ DÉCLARATION DES CONSTANTES __________________

a=3 #alpha
b=1/40 #beta
c=1/10 #gamma
d=2 #delta
e=1/20 #epsilon
f=1/40 #zeta
g=1/40 #eta
k=1/20 #theta
j=1/20 #iota
h=1 #kappa


# _________________ DÉCLARATION DES CRITÈRES D'ARRÊTS __________________
dt = 0.01
eps=1e-7
Nmax=100

# _________________ DÉCLARATION DE LA CONDITION INITIALE __________________

w=np.array([12,10,8])

# _________________ LANCEMENT DE LA SIMULATION __________________

temp=[]
xSol=[]
ySol=[]
zSol=[]
sol=w
# sol=u

t=0
tf=50
while (t<=tf):
    temp.append(t)
    (sol)=NewtonNDim(LoktaVolterra3,Jacobien3,sol,eps,Nmax,a,b,c,d,e,f,g,k,j,h,dt)
    xSol.append(sol[0])
    ySol.append(sol[1])
    zSol.append(sol[2])
    t=t+dt

plt.plot(temp,xSol,label="Lapins")
plt.plot(temp,ySol,label="Moutons")
plt.plot(temp,zSol,label="Loups")
plt.legend()
plt.xlabel("Temps")
plt.ylabel("Effectifs des populations en nombres d'individus")
plt.show()

