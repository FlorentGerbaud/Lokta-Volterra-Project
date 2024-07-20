from copy import deepcopy
import numpy as np
import time
import matplotlib.pyplot as plt
import math 



def LoktaVolterraRk2(z0,a,b,c,d):
    return np.array([a*z0[0]-b*z0[0]*z0[1],
                     c*z0[0]*z0[1]-d*z0[1]])
    
def LoktaVolterraRk3(z0,a,b,c,d,e,f,g,h,i,j):
    return np.array([z0[0]*(a-b*z0[0]-c*z0[1]),
                    z0[1]*(d-e*z0[1]-f*z0[0]-g*z0[2]),
                    z0[2]*(h*z0[1]-i*z0[2]-j)])
        

def RungeKutaR4_2(z0,h,a,b,c,d):
    k1=LoktaVolterraRk2(z0,a,b,c,d)
    cte=z0+(h/2)*k1
    k2=LoktaVolterraRk2(cte,a,b,c,d)
    cte2=z0+(h/2)*k2
    k3=LoktaVolterraRk2(cte2,a,b,c,d)
    cte3=z0+h*k3
    k4=LoktaVolterraRk2(cte3,a,b,c,d)
    # print("k2:",k2)
    return np.array([z0[0]+h*((k1[0]/6)+(k2[0]/3)+(k3[0]/3)+(k4[0]/6)),
                     z0[1]+h*((k1[1]/6)+(k2[1]/3)+(k3[1]/3)+(k4[1]/6))])
    
    
def RungeKutaR4_3(z0,a,b,c,d,e,f,g,h,i,j,dt):
    k1=LoktaVolterraRk3(z0,a,b,c,d,e,f,g,h,i,j)
    cte=z0+(dt/2)*k1
    k2=LoktaVolterraRk3(cte,a,b,c,d,e,f,g,h,i,j)
    cte2=z0+(dt/2)*k2
    k3=LoktaVolterraRk3(cte2,a,b,c,d,e,f,g,h,i,j)
    cte3=z0+dt*k3
    k4=LoktaVolterraRk3(cte3,a,b,c,d,e,f,g,h,i,j)
    # print("k2:",k2)
    return np.array([z0[0]+dt*((k1[0]/6)+(k2[0]/3)+(k3[0]/3)+(k4[0]/6)),
                     z0[1]+dt*((k1[1]/6)+(k2[1]/3)+(k3[1]/3)+(k4[1]/6)),
                     z0[2]+dt*((k1[2]/6)+(k2[2]/3)+(k3[2]/3)+(k4[2]/6))])

#_____________________ Parametre pour le système a 2 especes _____________________________#


#_____________________ Parametre pour le système a 3 especes _____________________________



temp=[]
xSol=[]
ySol=[]
zSol=[]
t=0
tf=2000


# ______________________________ CHOIX ESPECES _____________________________
choixSysEspece=2 #si 2 alors LV2 si 3 LV3
# __________________________________________________________________________


ens_x = []
ens_y = []

if(choixSysEspece==2):

    a = 2
    b = 1
    c = 0.2
    d = 2
    h = 0.1
    u = np.array([6,2])
    sol=u
    while (t<=tf):
        temp.append(t)
        (sol)=RungeKutaR4_2(sol,h,a,b,c,d)
        xSol.append(sol[0])
        ySol.append(sol[1])
        t=t+h

    plt.plot(temp,xSol,label="Lapins")
    plt.plot(temp,ySol,label="Loups")
    plt.grid()
    plt.legend()
    plt.xlabel("Temps")
    plt.ylabel("Effectif des populations en nombres d'individus")
    plt.show()

elif(choixSysEspece==3):
    a=3 # indice de l'espece x a determiner
    b=1/40
    c=1/10
    d=2 #indice de reproduction des prédateur
    e=1/20
    f=1/40
    g=1/20
    k=1/20 # indice reproduction des proies
    j=1/20
    h=1
    dt = 0.01
    w=np.array([6.,10.,4.])
    sol=w
    while (t<=tf):
        temp.append(t)
        (sol)=RungeKutaR4_3(sol,a,b,c,d,e,f,g,k,j,h,dt)
        # (sol)=RungeKutaR2_2(sol,dt,a,b,c,d)
        xSol.append(sol[0])
        ySol.append(sol[1])
        zSol.append(sol[2])
        t=t+dt

    plt.plot(temp,xSol,label="Lapins")
    plt.plot(temp,ySol,label="Moutons")
    plt.plot(temp,zSol,label="Loups")
    plt.grid()
    plt.legend()
    plt.xlabel("Temps")
    plt.ylabel("Effectif des populations en nombres d'individus")
    plt.show()
else:
    exit 



