import numpy as np
import matplotlib.pyplot as plt
from Riemann_modified import Riemann

c_max = 0.318
dx = 0.01
# case 1 - Sod problem
rhoL=1.0 
uL=0.0 
pL=1.0 
rhoR=0.125 
uR=0.0 
pR=0.1 
t_final = 0.25
L = 1.0
gamma = 1.4
ix = int(L/dx+1)
half = int(np.floor(ix/2))
x = np.linspace(0,L,ix)

# storage vectors
rho = np.zeros(ix)
u = np.zeros(ix)
p = np.zeros(ix)
e = np.zeros(ix)
Q = np.zeros((3,ix))
F_plus = np.zeros((3,ix))
F_minus = np.zeros((3,ix))
F_half = np.zeros((3,ix))

#establish initial condition vectors
rho[0:half] = rhoL
rho[half:ix] = rhoR
u[0:half] = uL
u[half:ix] = uR
p[0:half] = pL
p[half:ix] = pR
e = p/((gamma-1)*rho)
a = np.sqrt((gamma*p/rho))
M = u/a

t = 0
dt = c_max*dx/np.max(np.abs(u)+a)

#calculate initial Q and F vectors
Q[0,:] = rho
Q[1,:] = rho*u
Q[2,:] = rho*(e+u**2/2)

rhohalfL = np.zeros(ix)
rhohalfR = np.zeros(ix)
uhalfL = np.zeros(ix)
uhalfR = np.zeros(ix)
phalfL = np.zeros(ix)
phalfR = np.zeros(ix)

rhohalfL[:] = rho
rhohalfR[:] = rho
uhalfL[:] = u
uhalfR[:] = u
phalfL[:] = p
phalfR[:] = p

def MUSCL(vm1,v,vp1,vp2):
    
    kappa = -1
    rmp12 = (v-vm1+1.e-30)/(vp1-v+1.e-30)
    rmp32 = (vp1-v+1.e-30)/(vp2-vp1+1.e-30)
    phi_rmp12 =  (np.abs(rmp12)+rmp12)/(1+np.abs(rmp12))
    phi_rmp32 = (np.abs(rmp32)+rmp32)/(1+np.abs(rmp32))
    phi_rpp12 = (1/rmp32)*phi_rmp32
    phi_rpm12 = (1/rmp12)*phi_rmp12
    del_pm12 = (v-vm1)*phi_rpm12/2
    del_mp12 = (vp1-v)*phi_rmp12/2
    del_mp32 = (vp2-vp1)*phi_rmp32/2
    del_pp12 = (vp1-v)*phi_rpp12
    del_L = (1-kappa)*del_pm12+(1+kappa)*del_mp12
    del_R = (1-kappa)*del_mp32+(1+kappa)*del_pp12
    vL12 = v+0.5*del_L 
    vR12 = vp1-0.5*del_R
    
    return (vL12,vR12)     
    
while t <= t_final:
    for i in range(1,ix-2): # Loop for calculating Fn  
        
        im=i-1 
        ip=i+1
        ip2=i+2
        
        rhohalfL[i],rhohalfR[i] = MUSCL(rho[im],rho[i],rho[ip],rho[ip2])
        uhalfL[i],uhalfR[i] = MUSCL(u[im],u[i],u[ip],u[ip2])  
        phalfL[i],phalfR[i] = MUSCL(p[im],p[i],p[ip],p[ip2])  

    ahalfL = np.sqrt(gamma*phalfL/rhohalfL)
    ahalfR = np.sqrt(gamma*phalfR/rhohalfR)
    MhalfL = uhalfL/ahalfL
    MhalfR = uhalfR/ahalfR
    
    for i in range(0,ix): # Loop for calculating Fn 
                im=i-1 
                ip=i+1
                    
                # Neumann BCs implemented here using a series of if statements
                if i == 0:
                    im = 1
                if i == ix-1:
                    ip = ix-2
                             
                F_plus[0,i] = ((1/4)*rhohalfL[i]*ahalfL[i]*(1+MhalfL[i])**2)
                F_plus[1,i] = ((1/4)*rhohalfL[i]*ahalfL[i]*(1+MhalfL[i])**2)*(2*ahalfL[i]/gamma)*(((gamma-1)/2)*MhalfL[i]+1)
                F_plus[2,i] = ((1/4)*rhohalfL[i]*ahalfL[i]*(1+MhalfL[i])**2)*(2*ahalfL[i]**2/(gamma**2-1))*(((gamma-1)/2)*MhalfL[i]+1)**2
                
                F_minus[0,i] = -((1/4)*rhohalfR[i]*ahalfR[i]*(1-MhalfR[i])**2)
                F_minus[1,i] = -((1/4)*rhohalfR[i]*ahalfR[i]*(1-MhalfR[i])**2)*(2*ahalfR[i]/gamma)*(((gamma-1)/2)*MhalfR[i]-1)
                F_minus[2,i] = -((1/4)*rhohalfR[i]*ahalfR[i]*(1-MhalfR[i])**2)*(2*ahalfR[i]**2/(gamma**2-1))*(((gamma-1)/2)*MhalfR[i]-1)**2
            
                F_half[:,i] = F_plus[:,i] + F_minus[:,i]
        
    for i in range(0,ix): # Q_update
        im=i-1 
        ip=i+1        
        # Neumann BCs implemented here using a series of if statements
        if i == 0:
            im = 1
        if i == ix-1:
            ip = ix-2
                
        Q[:,i] = Q[:,i]-(dt/dx)*(F_half[:,i]-F_half[:,im])
    
    rho = Q[0,:]
    u = Q[1,:]/rho
    e = Q[2,:]/rho-u**2/2
    p = e*(gamma-1)*rho
    a = np.sqrt((gamma*p/rho))
    M = u/a
    
    dt = c_max*dx/np.max(np.abs(u)+a)
    t = t+dt
    
    print(t)

if __name__ == '__main__':
 R=Riemann()
 rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
 case=1 # case numbers are summarized above 
 gam=1.4  
 [rhoexact,uexact,pexact,eexact,xexact] = R.exact(gam,max_time[case],case,"rusanov"+str(case)) 


fig1, axs1 = plt.subplots(4,figsize=(6,8))
fig1.suptitle("Comparison of shock-capturing methods, Test 1, $Δx = 0.025$, $c_{max} = 0.5$",fontsize=10)

axs1[0].plot(xexact, rhoexact,"--k")
axs1[1].plot(xexact, uexact,"--k")
axs1[2].plot(xexact, pexact,"--k")
axs1[3].plot(xexact, eexact,"--k",label="Analytical")

axs1[0].plot(x, rho)
axs1[1].plot(x, u)
axs1[2].plot(x, p)
axs1[3].plot(x, e,label="Van Leer MUSCL")