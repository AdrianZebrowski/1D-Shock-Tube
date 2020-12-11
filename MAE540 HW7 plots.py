#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 01:58:56 2020

@author: adrianzebrowski
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from Riemann_modified import Riemann
from vanleer1D import vanleer1D
from roe1D import roe1D
from vanleer1DMUSCL import vanleer1DMUSCL
from roe1DMUSCL import roe1DMUSCL

if __name__ == '__main__':
 R=Riemann()
 rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
 case=1 # case numbers are summarized above 
 gam=1.4  
 [rhoexact,uexact,pexact,eexact,xexact] = R.exact(gam,max_time[case],case,"rusanov"+str(case)) 

### PART 1
fig1, axs1 = plt.subplots(4,figsize=(6,8))
fig1.suptitle("Comparison of shock-capturing methods, Test 1, $Δx = 0.025$, $c_{max} = 0.5$",fontsize=10)

[rho1,u1,p1,e1,x1] = vanleer1D(1,0.5,0.025)
[rho2,u2,p2,e2,x2] = roe1D(1,0.5,0.025)

axs1[0].plot(xexact, rhoexact,"--k")
axs1[1].plot(xexact, uexact,"--k")
axs1[2].plot(xexact, pexact,"--k")
axs1[3].plot(xexact, eexact,"--k",label="Analytical")

axs1[0].plot(x1, rho1)
axs1[0].grid()
axs1[1].plot(x1, u1)
axs1[1].grid()
axs1[2].plot(x1, p1)
axs1[2].grid()
axs1[3].plot(x1, e1,label="Van leer")
axs1[3].grid()

axs1[0].plot(x2, rho2)
axs1[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs1[0].set_xlim([0, 1])
axs1[1].plot(x2, u2)
axs1[1].set(ylabel="$u$ $(m/s)$")
axs1[1].set_xlim([0, 1])
axs1[1].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs1[2].plot(x2, p2)
axs1[2].set(ylabel="$p$ $(Pa)$")
axs1[2].set_xlim([0, 1])
axs1[2].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs1[3].plot(x2, e2,label="Roe")
axs1[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs1[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs1[3].set_xlim([0, 1])
axs1[3].yaxis.set_major_locator(ticker.MultipleLocator(0.2))

fig1.subplots_adjust(top=0.95)
fig1.savefig('HW7 Part 1.png', dpi=300)

##########

### PART 2
fig2, axs1 = plt.subplots(4,figsize=(6,8))
fig2.suptitle("Comparison of shock-capturing methods, Test 1, $Δx = 0.0125$, $c_{max} = 0.5$",fontsize=10)

[rho1,u1,p1,e1,x1] = vanleer1D(1,0.5,0.0125)
[rho2,u2,p2,e2,x2] = roe1D(1,0.5,0.0125)

axs1[0].plot(xexact, rhoexact,"--k")
axs1[1].plot(xexact, uexact,"--k")
axs1[2].plot(xexact, pexact,"--k")
axs1[3].plot(xexact, eexact,"--k",label="Analytical")

axs1[0].plot(x1, rho1)
axs1[0].grid()
axs1[1].plot(x1, u1)
axs1[1].grid()
axs1[2].plot(x1, p1)
axs1[2].grid()
axs1[3].plot(x1, e1,label="Van leer")
axs1[3].grid()

axs1[0].plot(x2, rho2)
axs1[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs1[0].set_xlim([0, 1])
axs1[1].plot(x2, u2)
axs1[1].set(ylabel="$u$ $(m/s)$")
axs1[1].set_xlim([0, 1])
axs1[1].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs1[2].plot(x2, p2)
axs1[2].set(ylabel="$p$ $(Pa)$")
axs1[2].set_xlim([0, 1])
axs1[2].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs1[3].plot(x2, e2,label="Roe")
axs1[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs1[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs1[3].set_xlim([0, 1])
axs1[3].yaxis.set_major_locator(ticker.MultipleLocator(0.2))

fig2.subplots_adjust(top=0.95)
fig2.savefig('HW7 Part 2.png', dpi=300)

###############

#### PART 3
dx_array = np.array([0.1,0.05,0.025,0.0125,0.00625,0.003125,0.0015625])

L2rho = np.zeros((2,len(dx_array)))
L2u = np.zeros((2,len(dx_array)))
L2p = np.zeros((2,len(dx_array)))
L2e = np.zeros((2,len(dx_array)))

first_order = np.zeros((4,len(dx_array)))
second_order = np.zeros((4,len(dx_array)))

for i, dx in enumerate(dx_array):     
    
    [rho1,u1,p1,e1,x1] = vanleer1D(1,0.5,dx)
    [rho2,u2,p2,e2,x2] = roe1D(1,0.5,dx)
    
    rhoexactinterp = np.interp(x1,xexact,rhoexact)
    uexactinterp = np.interp(x1,xexact,uexact)
    pexactinterp = np.interp(x1,xexact,pexact)
    eexactinterp = np.interp(x1,xexact,eexact)

    L2rho[0,i] = np.linalg.norm(rho1-rhoexactinterp)/np.linalg.norm(rhoexactinterp)
    L2rho[1,i] = np.linalg.norm(rho2-rhoexactinterp)/np.linalg.norm(rhoexactinterp)
    
    L2u[0,i] = np.linalg.norm(u1-uexactinterp)/np.linalg.norm(uexactinterp) 
    L2u[1,i] = np.linalg.norm(u2-uexactinterp)/np.linalg.norm(uexactinterp)  
    
    L2p[0,i] = np.linalg.norm(p1-pexactinterp)/np.linalg.norm(pexactinterp) 
    L2p[1,i] = np.linalg.norm(p2-pexactinterp)/np.linalg.norm(pexactinterp) 
    
    L2e[0,i] = np.linalg.norm(e1-eexactinterp)/np.linalg.norm(eexactinterp) 
    L2e[1,i] = np.linalg.norm(e2-eexactinterp)/np.linalg.norm(eexactinterp)  
    
    first_order[0,i] = L2rho[1,0]/(2**i)
    second_order[0,i] = L2rho[1,0]/(2**(2*i))
    first_order[1,i] = L2u[1,0]/(2**i)
    second_order[1,i] = L2u[1,0]/(2**(2*i))
    first_order[2,i] = L2p[1,0]/(2**i)
    second_order[2,i] = L2p[1,0]/(2**(2*i))
    first_order[3,i] = L2e[1,0]/(2**i)
    second_order[3,i] = L2e[1,0]/(2**(2*i))

fig3, axs4 = plt.subplots(4,figsize=(6,8))
fig3.suptitle("Comparison of L2 error norms, Test 1, $c_{max} = 0.5$",fontsize=10)

axs4[0].loglog(1/dx_array,first_order[0,:],"k")
axs4[1].loglog(1/dx_array,first_order[1,:],"k")
axs4[2].loglog(1/dx_array,first_order[2,:],"k")
axs4[3].loglog(1/dx_array,first_order[3,:],"k",label="$Δx$ convergence rate")

axs4[0].loglog(1/dx_array,second_order[0,:],"--k")
axs4[1].loglog(1/dx_array,second_order[1,:],"--k")
axs4[2].loglog(1/dx_array,second_order[2,:],"--k")
axs4[3].loglog(1/dx_array,second_order[3,:],"--k",label="$(Δx)^2$ convergence rate")

axs4[0].loglog(1/dx_array,L2rho[0,:],"o")
axs4[0].grid()
axs4[1].loglog(1/dx_array,L2u[0,:],"o")
axs4[1].grid()
axs4[2].loglog(1/dx_array,L2p[0,:],"o")
axs4[2].grid()
axs4[3].loglog(1/dx_array,L2e[0,:],"o",label="Van Leer")
axs4[3].grid()


axs4[0].loglog(1/dx_array,L2rho[1,:],"^")
axs4[0].set(ylabel="$log(\epsilon)$, $ρ$")
axs4[0].set_xlim([9, 1e3])
axs4[1].loglog(1/dx_array,L2u[1,:],"^")
axs4[1].set(ylabel="$log(\epsilon)$, $u$")
axs4[1].set_xlim([9, 1e3])
axs4[2].loglog(1/dx_array,L2p[1,:],"^")
axs4[2].set(ylabel="$log(\epsilon)$, $p$")
axs4[2].set_xlim([9, 1e3])
axs4[3].loglog(1/dx_array,L2e[1,:],"^",label="Roe")
axs4[3].set(xlabel="$log(1/Δx)$",ylabel="$log(\epsilon)$, $e$")
axs4[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.35),ncol=3)
axs4[3].set_xlim([9, 1e3])
fig3.subplots_adjust(top=0.95)
fig3.savefig('HW7 Part 3.png', dpi=300)

#### Part 5
### 123 problem

if __name__ == '__main__':
 R=Riemann()
 rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
 case=2 # case numbers are summarized above 
 gam=1.4  
 [rhoexact,uexact,pexact,eexact,xexact] = R.exact(gam,max_time[case],case,"rusanov"+str(case)) 

fig4, axs5 = plt.subplots(4,figsize=(6,8))
fig4.suptitle("Comparison of shock-capturing methods, Test 2, $Δx = 0.0125$, $c_{max} = 0.5$",fontsize=10)

[rho1,u1,p1,e1,x1] = vanleer1D(2,0.5,0.0125)
[rho2,u2,p2,e2,x2] = roe1D(2,0.5,0.0125)

axs5[0].plot(xexact, rhoexact,"--k")
axs5[1].plot(xexact, uexact,"--k")
axs5[2].plot(xexact, pexact,"--k")
axs5[3].plot(xexact, eexact,"--k",label="Analytical")

axs5[0].plot(x1, rho1)
axs5[0].grid()
axs5[1].plot(x1, u1)
axs5[1].grid()
axs5[2].plot(x1, p1)
axs5[2].grid()
axs5[3].plot(x1, e1,label="Van Leer")
axs5[3].grid()

axs5[0].plot(x2, rho2)
axs5[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs5[0].set_xlim([0, 1])
axs5[1].plot(x2, u2)
axs5[1].set(ylabel="$u$ $(m/s)$")
axs5[1].set_xlim([0, 1])
axs5[2].plot(x2, p2)
axs5[2].set(ylabel="$p$ $(Pa)$")
axs5[2].set_xlim([0, 1])
axs5[3].plot(x2, e2,label="Roe")
axs5[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs5[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs5[3].set_xlim([0, 1])

fig4.subplots_adjust(top=0.95)
fig4.savefig('HW7 Part 4_123.png', dpi=300)

##########
if __name__ == '__main__':
 R=Riemann()
 rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
 case=3 # case numbers are summarized above 
 gam=1.4  
 [rhoexact,uexact,pexact,eexact,xexact] = R.exact(gam,max_time[case],case,"rusanov"+str(case)) 

fig5, axs5 = plt.subplots(4,figsize=(6,8))
fig5.suptitle("Comparison of shock-capturing methods, Test 3, $Δx = 0.0125$, $c_{max} = 0.5$",fontsize=10)

[rho1,u1,p1,e1,x1] = vanleer1D(3,0.5,0.0125)
[rho2,u2,p2,e2,x2] = roe1D(3,0.5,0.0125)

axs5[0].plot(xexact, rhoexact,"--k")
axs5[1].plot(xexact, uexact,"--k")
axs5[2].plot(xexact, pexact,"--k")
axs5[3].plot(xexact, eexact,"--k",label="Analytical")

axs5[0].plot(x1, rho1)
axs5[0].grid()
axs5[1].plot(x1, u1)
axs5[1].grid()
axs5[2].plot(x1, p1)
axs5[2].grid()
axs5[3].plot(x1, e1,label="Van Leer")
axs5[3].grid()

axs5[0].plot(x2, rho2)
axs5[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs5[0].set_xlim([0, 1])
axs5[1].plot(x2, u2)
axs5[1].set(ylabel="$u$ $(m/s)$")
axs5[1].set_xlim([0, 1])
axs5[2].plot(x2, p2)
axs5[2].set(ylabel="$p$ $(Pa)$")
axs5[2].set_xlim([0, 1])
axs5[3].plot(x2, e2,label="Roe")
axs5[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs5[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs5[3].set_xlim([0, 1])

fig5.subplots_adjust(top=0.95)
fig5.savefig('HW7 Part 4_blast1.png', dpi=300)

####
if __name__ == '__main__':
 R=Riemann()
 rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
 case=4 # case numbers are summarized above 
 gam=1.4  
 [rhoexact,uexact,pexact,eexact,xexact] = R.exact(gam,max_time[case],case,"rusanov"+str(case)) 

fig6, axs5 = plt.subplots(4,figsize=(6,8))
fig6.suptitle("Comparison of shock-capturing methods, Test 4, $Δx = 0.0125$, $c_{max} = 0.5$",fontsize=10)

[rho1,u1,p1,e1,x1] = vanleer1D(4,0.5,0.0125)
[rho2,u2,p2,e2,x2] = roe1D(4,0.5,0.0125)

axs5[0].plot(xexact, rhoexact,"--k")
axs5[1].plot(xexact, uexact,"--k")
axs5[2].plot(xexact, pexact,"--k")
axs5[3].plot(xexact, eexact,"--k",label="Analytical")

axs5[0].plot(x1, rho1)
axs5[0].grid()
axs5[1].plot(x1, u1)
axs5[1].grid()
axs5[2].plot(x1, p1)
axs5[2].grid()
axs5[3].plot(x1, e1,label="Van Leer")
axs5[3].grid()

axs5[0].plot(x2, rho2)
axs5[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs5[0].set_xlim([0, 1])
axs5[1].plot(x2, u2)
axs5[1].set(ylabel="$u$ $(m/s)$")
axs5[1].set_xlim([0, 1])
axs5[2].plot(x2, p2)
axs5[2].set(ylabel="$p$ $(Pa)$")
axs5[2].set_xlim([0, 1])
axs5[3].plot(x2, e2,label="Roe")
axs5[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs5[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs5[3].set_xlim([0, 1])

fig6.subplots_adjust(top=0.95)
fig6.savefig('HW7 Part 4_blast2.png', dpi=300)

#####
if __name__ == '__main__':
 R=Riemann()
 rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
 case=5 # case numbers are summarized above 
 gam=1.4  
 [rhoexact,uexact,pexact,eexact,xexact] = R.exact(gam,max_time[case],case,"rusanov"+str(case)) 

fig7, axs5 = plt.subplots(4,figsize=(6,8))
fig7.suptitle("Comparison of shock-capturing methods, Test 5, $Δx = 0.0125$, $c_{max} = 0.5$",fontsize=10)

[rho1,u1,p1,e1,x1] = vanleer1D(5,0.5,0.0125)
[rho2,u2,p2,e2,x2] = roe1D(5,0.5,0.0125)

axs5[0].plot(xexact, rhoexact,"--k")
axs5[1].plot(xexact, uexact,"--k")
axs5[2].plot(xexact, pexact,"--k")
axs5[3].plot(xexact, eexact,"--k",label="Analytical")

axs5[0].plot(x1, rho1)
axs5[0].grid()
axs5[1].plot(x1, u1)
axs5[1].grid()
axs5[2].plot(x1, p1)
axs5[2].grid()
axs5[3].plot(x1, e1,label="Van Leer")
axs5[3].grid()

axs5[0].plot(x2, rho2)
axs5[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs5[0].set_xlim([0, 1])
axs5[1].plot(x2, u2)
axs5[1].set(ylabel="$u$ $(m/s)$")
axs5[1].set_xlim([0, 1])
axs5[2].plot(x2, p2)
axs5[2].set(ylabel="$p$ $(Pa)$")
axs5[2].set_xlim([0, 1])
axs5[3].plot(x2, e2,label="Roe")
axs5[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs5[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs5[3].set_xlim([0, 1])

fig7.subplots_adjust(top=0.95)
fig7.savefig('HW7 Part 4_shockcollision.png', dpi=300)