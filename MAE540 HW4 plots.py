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
from functions import maccormackFD, lax_fredrich, rusanov, maccormackFV, maccormackFCT

if __name__ == '__main__':
 R=Riemann()
 rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
 case=1 # case numbers are summarized above 
 gam=1.4  
 [rhoexact,uexact,pexact,eexact,xexact] = R.exact(gam,max_time[case],case,"rusanov"+str(case)) 

### PART 1
fig1, axs1 = plt.subplots(4,figsize=(6,8))
fig1.suptitle("Comparison of shock-capturing methods, Test 1, $Δx = 0.025$, $c_{max} = 1.0$",fontsize=10)

[rho1,u1,p1,e1,x1] = maccormackFD("Case 1",1.0,0.025,1.4)
[rho2,u2,p2,e2,x2] = lax_fredrich("Case 1",1.0,0.025,1.4)
[rho3,u3,p3,e3,x3] = rusanov("Case 1",1.0,0.025,1.4)

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
axs1[3].plot(x1, e1,label="MacCormack FD")
axs1[3].grid()

axs1[0].plot(x2, rho2)
axs1[1].plot(x2, u2)
axs1[2].plot(x2, p2)
axs1[3].plot(x2, e2,label="Lax-Fredrich")

axs1[0].plot(x3, rho3)
axs1[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs1[0].set_xlim([0, 1])
axs1[1].plot(x3, u3)
axs1[1].set(ylabel="$u$ $(m/s)$")
axs1[1].set_xlim([0, 1])
axs1[1].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs1[2].plot(x3, p3)
axs1[2].set(ylabel="$p$ $(Pa)$")
axs1[2].set_xlim([0, 1])
axs1[2].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs1[3].plot(x3, e3,label="Rusanov")
axs1[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs1[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs1[3].set_xlim([0, 1])
axs1[3].yaxis.set_major_locator(ticker.MultipleLocator(0.2))

fig1.subplots_adjust(top=0.95)
fig1.savefig('HW4 Part 1.png', dpi=300)

##########

### PART 2
fig2, axs2 = plt.subplots(4,figsize=(6,8))
fig2.suptitle("Comparison of shock-capturing methods, Test 1, $Δx = 0.025$, $c_{max} = 0.8$",fontsize=10)

[rho1,u1,p1,e1,x1] = maccormackFD("Case 1",0.8,0.025,1.4)
[rho2,u2,p2,e2,x2] = lax_fredrich("Case 1",0.8,0.025,1.4)
[rho3,u3,p3,e3,x3] = rusanov("Case 1",0.8,0.025,1.4)

axs2[0].plot(xexact, rhoexact,"--k")
axs2[1].plot(xexact, uexact,"--k")
axs2[2].plot(xexact, pexact,"--k")
axs2[3].plot(xexact, eexact,"--k",label="Analytical")

axs2[0].plot(x1, rho1)
axs2[0].grid()
axs2[1].plot(x1, u1)
axs2[1].grid()
axs2[2].plot(x1, p1)
axs2[2].grid()
axs2[3].plot(x1, e1,label="MacCormack FD")
axs2[3].grid()

axs2[0].plot(x2, rho2)
axs2[1].plot(x2, u2)
axs2[2].plot(x2, p2)
axs2[3].plot(x2, e2,label="Lax-Fredrich")

axs2[0].plot(x3, rho3)
axs2[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs2[0].set_xlim([0, 1])
axs2[1].plot(x3, u3)
axs2[1].set(ylabel="$u$ $(m/s)$")
axs2[1].set_xlim([0, 1])
axs2[1].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs2[2].plot(x3, p3)
axs2[2].set(ylabel="$p$ $(Pa)$")
axs2[2].set_xlim([0, 1])
axs2[2].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs2[3].plot(x3, e3,label="Rusanov")
axs2[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs2[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs2[3].set_xlim([0, 1])
axs2[3].yaxis.set_major_locator(ticker.MultipleLocator(0.2))

fig2.subplots_adjust(top=0.95)
fig2.savefig('HW4 Part 2.png', dpi=300)

###############

### PART 3
fig3, axs3 = plt.subplots(4,figsize=(6,8))
fig3.suptitle("Comparison of shock-capturing methods, Test 1, $Δx = 0.0125$, $c_{max} = 1.0$",fontsize=10)

[rho1,u1,p1,e1,x1] = maccormackFD("Case 1",1.0,0.0125,1.4)
[rho2,u2,p2,e2,x2] = lax_fredrich("Case 1",1.0,0.0125,1.4)
[rho3,u3,p3,e3,x3] = rusanov("Case 1",1.0,0.0125,1.4)

axs3[0].plot(xexact, rhoexact,"--k")
axs3[1].plot(xexact, uexact,"--k")
axs3[2].plot(xexact, pexact,"--k")
axs3[3].plot(xexact, eexact,"--k",label="Analytical")

axs3[0].plot(x1, rho1)
axs3[0].grid()
axs3[1].plot(x1, u1)
axs3[1].grid()
axs3[2].plot(x1, p1)
axs3[2].grid()
axs3[3].plot(x1, e1,label="MacCormack FD")
axs3[3].grid()

axs3[0].plot(x2, rho2)
axs3[1].plot(x2, u2)
axs3[2].plot(x2, p2)
axs3[3].plot(x2, e2,label="Lax-Fredrich")

axs3[0].plot(x3, rho3)
axs3[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs3[0].set_xlim([0, 1])
axs3[1].plot(x3, u3)
axs3[1].set(ylabel="$u$ $(m/s)$")
axs3[1].set_xlim([0, 1])
axs3[1].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs3[2].plot(x3, p3)
axs3[2].set(ylabel="$p$ $(Pa)$")
axs3[2].set_xlim([0, 1])
axs3[2].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs3[3].plot(x3, e3,label="Rusanov")
axs3[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs3[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs3[3].set_xlim([0, 1])
axs3[3].yaxis.set_major_locator(ticker.MultipleLocator(0.2))

fig3.subplots_adjust(top=0.95)
fig3.savefig('HW4 Part 3.png', dpi=300)

####

#### PART 4
dx_array = np.array([0.1,0.05,0.025,0.0125,0.00625,0.003125,0.0015625])

L2rho = np.zeros((3,len(dx_array)))
L2u = np.zeros((3,len(dx_array)))
L2p = np.zeros((3,len(dx_array)))
L2e = np.zeros((3,len(dx_array)))

first_order = np.zeros((4,len(dx_array)))
second_order = np.zeros((4,len(dx_array)))

for i, dx in enumerate(dx_array):     
    
    [rho1,u1,p1,e1,x1] = maccormackFD("Case 1",1.0,dx,1.4)
    [rho2,u2,p2,e2,x2] = lax_fredrich("Case 1",1.0,dx,1.4)
    [rho3,u3,p3,e3,x3] = rusanov("Case 1",1.0,dx,1.4)
    
    rhoexactinterp = np.interp(x1,xexact,rhoexact)
    uexactinterp = np.interp(x1,xexact,uexact)
    pexactinterp = np.interp(x1,xexact,pexact)
    eexactinterp = np.interp(x1,xexact,eexact)

    L2rho[0,i] = np.linalg.norm(rho1-rhoexactinterp)/np.linalg.norm(rhoexactinterp)
    L2rho[1,i] = np.linalg.norm(rho2-rhoexactinterp)/np.linalg.norm(rhoexactinterp)
    L2rho[2,i] = np.linalg.norm(rho3-rhoexactinterp)/np.linalg.norm(rhoexactinterp)
    
    L2u[0,i] = np.linalg.norm(u1-uexactinterp)/np.linalg.norm(uexactinterp) 
    L2u[1,i] = np.linalg.norm(u2-uexactinterp)/np.linalg.norm(uexactinterp) 
    L2u[2,i] = np.linalg.norm(u3-uexactinterp)/np.linalg.norm(uexactinterp)  
    
    L2p[0,i] = np.linalg.norm(p1-pexactinterp)/np.linalg.norm(pexactinterp) 
    L2p[1,i] = np.linalg.norm(p2-pexactinterp)/np.linalg.norm(pexactinterp) 
    L2p[2,i] = np.linalg.norm(p3-pexactinterp)/np.linalg.norm(pexactinterp) 
    
    L2e[0,i] = np.linalg.norm(e1-eexactinterp)/np.linalg.norm(eexactinterp) 
    L2e[1,i] = np.linalg.norm(e2-eexactinterp)/np.linalg.norm(eexactinterp)  
    L2e[2,i] = np.linalg.norm(e3-eexactinterp)/np.linalg.norm(eexactinterp)  
    
    first_order[0,i] = L2rho[1,0]/(2**i)
    second_order[0,i] = L2rho[1,0]/(2**(2*i))
    first_order[1,i] = L2u[1,0]/(2**i)
    second_order[1,i] = L2u[1,0]/(2**(2*i))
    first_order[2,i] = L2p[1,0]/(2**i)
    second_order[2,i] = L2p[1,0]/(2**(2*i))
    first_order[3,i] = L2e[1,0]/(2**i)
    second_order[3,i] = L2e[1,0]/(2**(2*i))

fig4, axs4 = plt.subplots(4,figsize=(6,8.5))
fig4.suptitle("Comparison of L2 error norms, Test 1, $c_{max} = 1.0$",fontsize=10)

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
axs4[3].loglog(1/dx_array,L2e[0,:],"o",label="MacCormack FD")
axs4[3].grid()

axs4[0].loglog(1/dx_array,L2rho[1,:],"s")
axs4[1].loglog(1/dx_array,L2u[1,:],"s")
axs4[2].loglog(1/dx_array,L2p[1,:],"s")
axs4[3].loglog(1/dx_array,L2e[1,:],"s",label="Lax-Fredrich")

axs4[0].loglog(1/dx_array,L2rho[2,:],"^")
axs4[0].set(ylabel="$log(\epsilon)$, $ρ$")
axs4[0].set_xlim([9, 1e3])
axs4[1].loglog(1/dx_array,L2u[2,:],"^")
axs4[1].set(ylabel="$log(\epsilon)$, $u$")
axs4[1].set_xlim([9, 1e3])
axs4[2].loglog(1/dx_array,L2p[2,:],"^")
axs4[2].set(ylabel="$log(\epsilon)$, $p$")
axs4[2].set_xlim([9, 1e3])
axs4[3].loglog(1/dx_array,L2e[2,:],"^",label="Rusanov")
axs4[3].set(xlabel="$log(1/Δx)$",ylabel="$log(\epsilon)$, $e$")
axs4[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.35),ncol=3)
axs4[3].set_xlim([9, 1e3])

fig4.subplots_adjust(top=0.95)
fig4.savefig('HW4 Part 4.png', dpi=300)

#### Part 5
### 123 problem

if __name__ == '__main__':
 R=Riemann()
 rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
 case=2 # case numbers are summarized above 
 gam=1.4  
 [rhoexact,uexact,pexact,eexact,xexact] = R.exact(gam,max_time[case],case,"rusanov"+str(case)) 

fig5, axs5 = plt.subplots(4,figsize=(6,8))
fig5.suptitle("Comparison of shock-capturing methods, Test 2, $Δx = 0.0125$, $c_{max} = 1.0$",fontsize=10)

[rho1,u1,p1,e1,x1] = maccormackFD("Case 2",1.0,0.0125,1.4)
[rho2,u2,p2,e2,x2] = lax_fredrich("Case 2",1.0,0.0125,1.4)
[rho3,u3,p3,e3,x3] = rusanov("Case 2",1.0,0.0125,1.4)

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
axs5[3].plot(x1, e1,label="MacCormack FD")
axs5[3].grid()

axs5[0].plot(x2, rho2)
axs5[1].plot(x2, u2)
axs5[2].plot(x2, p2)
axs5[3].plot(x2, e2,label="Lax-Fredrich")

axs5[0].plot(x3, rho3)
axs5[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs5[0].set_xlim([0, 1])
axs5[1].plot(x3, u3)
axs5[1].set(ylabel="$u$ $(m/s)$")
axs5[1].set_xlim([0, 1])
axs5[2].plot(x3, p3)
axs5[2].set(ylabel="$p$ $(Pa)$")
axs5[2].set_xlim([0, 1])
axs5[3].plot(x3, e3,label="Rusanov")
axs5[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs5[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs5[3].set_xlim([0, 1])

fig5.subplots_adjust(top=0.95)
fig5.savefig('HW4 Part 5_123.png', dpi=300)

##########
if __name__ == '__main__':
 R=Riemann()
 rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
 case=3 # case numbers are summarized above 
 gam=1.4  
 [rhoexact,uexact,pexact,eexact,xexact] = R.exact(gam,max_time[case],case,"rusanov"+str(case)) 
 
fig6, axs6 = plt.subplots(4,figsize=(6,8))
fig6.suptitle("Comparison of shock-capturing methods, Test 3, $Δx = 0.0125$, $c_{max} = 1.0$",fontsize=10)

[rho1,u1,p1,e1,x1] = maccormackFD("Case 3",1.0,0.0125,1.4)
[rho2,u2,p2,e2,x2] = lax_fredrich("Case 3",1.0,0.0125,1.4)
[rho3,u3,p3,e3,x3] = rusanov("Case 3",1.0,0.0125,1.4)

axs6[0].plot(xexact, rhoexact,"--k")
axs6[1].plot(xexact, uexact,"--k")
axs6[2].plot(xexact, pexact,"--k")
axs6[3].plot(xexact, eexact,"--k",label="Analytical")

axs6[0].plot(x1, rho1)
axs6[0].grid()
axs6[1].plot(x1, u1)
axs6[1].grid()
axs6[2].plot(x1, p1)
axs6[2].grid()
axs6[3].plot(x1, e1,label="MacCormack FD")
axs6[3].grid()

axs6[0].plot(x2, rho2)
axs6[1].plot(x2, u2)
axs6[2].plot(x2, p2)
axs6[3].plot(x2, e2,label="Lax-Fredrich")

axs6[0].plot(x3, rho3)
axs6[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs6[0].set_xlim([0, 1])
axs6[1].plot(x3, u3)
axs6[1].set(ylabel="$u$ $(m/s)$")
axs6[1].set_xlim([0, 1])
axs6[2].plot(x3, p3)
axs6[2].set(ylabel="$p$ $(Pa)$")
axs6[2].set_xlim([0, 1])
axs6[3].plot(x3, e3,label="Rusanov")
axs6[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs6[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs6[3].set_xlim([0, 1])

fig6.subplots_adjust(top=0.95)
fig6.savefig('HW4 Part 5_blast1.png', dpi=300)

####
if __name__ == '__main__':
 R=Riemann()
 rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
 case=4 # case numbers are summarized above 
 gam=1.4  
 [rhoexact,uexact,pexact,eexact,xexact] = R.exact(gam,max_time[case],case,"rusanov"+str(case)) 
 
fig7, axs7 = plt.subplots(4,figsize=(6,8))
fig7.suptitle("Comparison of shock-capturing methods, Test 4, $Δx = 0.0125$, $c_{max} = 1.0$",fontsize=10)

[rho1,u1,p1,e1,x1] = maccormackFD("Case 4",1.0,0.0125,1.4)
[rho2,u2,p2,e2,x2] = lax_fredrich("Case 4",1.0,0.0125,1.4)
[rho3,u3,p3,e3,x3] = rusanov("Case 4",1.0,0.0125,1.4)

axs7[0].plot(xexact, rhoexact,"--k")
axs7[1].plot(xexact, uexact,"--k")
axs7[2].plot(xexact, pexact,"--k")
axs7[3].plot(xexact, eexact,"--k",label="Analytical")

axs7[0].plot(x1, rho1)
axs7[0].grid()
axs7[1].plot(x1, u1)
axs7[1].grid()
axs7[2].plot(x1, p1)
axs7[2].grid()
axs7[3].plot(x1, e1,label="MacCormack FD")
axs7[3].grid()

axs7[0].plot(x2, rho2)
axs7[1].plot(x2, u2)
axs7[2].plot(x2, p2)
axs7[3].plot(x2, e2,label="Lax-Fredrich")

axs7[0].plot(x3, rho3)
axs7[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs7[0].set_xlim([0, 1])
axs7[1].plot(x3, u3)
axs7[1].set(ylabel="$u$ $(m/s)$")
axs7[1].set_xlim([0, 1])
axs7[2].plot(x3, p3)
axs7[2].set(ylabel="$p$ $(Pa)$")
axs7[2].set_xlim([0, 1])
axs7[3].plot(x3, e3,label="Rusanov")
axs7[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs7[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs7[3].set_xlim([0, 1])

fig7.subplots_adjust(top=0.95)
fig7.savefig('HW4 Part 5_blast2.png', dpi=300)

#####
if __name__ == '__main__':
 R=Riemann()
 rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
 case=5 # case numbers are summarized above 
 gam=1.4  
 [rhoexact,uexact,pexact,eexact,xexact] = R.exact(gam,max_time[case],case,"rusanov"+str(case)) 
 
fig8, axs8 = plt.subplots(4,figsize=(6,8))
fig8.suptitle("Comparison of shock-capturing methods, Test 5, $Δx = 0.0125$, $c_{max} = 1.0$",fontsize=10)

[rho1,u1,p1,e1,x1] = maccormackFD("Case 5",1.0,0.0125,1.4)
[rho2,u2,p2,e2,x2] = lax_fredrich("Case 5",1.0,0.0125,1.4)
[rho3,u3,p3,e3,x3] = rusanov("Case 5",1.0,0.0125,1.4)

axs8[0].plot(xexact, rhoexact,"--k")
axs8[1].plot(xexact, uexact,"--k")
axs8[2].plot(xexact, pexact,"--k")
axs8[3].plot(xexact, eexact,"--k",label="Analytical")

axs8[0].plot(x1, rho1)
axs8[0].grid()
axs8[1].plot(x1, u1)
axs8[1].grid()
axs8[2].plot(x1, p1)
axs8[2].grid()
axs8[3].plot(x1, e1,label="MacCormack FD")
axs8[3].grid()

axs8[0].plot(x2, rho2)
axs8[1].plot(x2, u2)
axs8[2].plot(x2, p2)
axs8[3].plot(x2, e2,label="Lax-Fredrich")

axs8[0].plot(x3, rho3)
axs8[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs8[0].set_xlim([0, 1])
axs8[1].plot(x3, u3)
axs8[1].set(ylabel="$u$ $(m/s)$")
axs8[1].set_xlim([0, 1])
axs8[2].plot(x3, p3)
axs8[2].set(ylabel="$p$ $(Pa)$")
axs8[2].set_xlim([0, 1])
axs8[3].plot(x3, e3,label="Rusanov")
axs8[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs8[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs8[3].set_xlim([0, 1])

fig8.subplots_adjust(top=0.95)
fig8.savefig('HW4 Part 5_shockcollision.png', dpi=300)

### PART 6.1
if __name__ == '__main__':
 R=Riemann()
 rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
 case=1 # case numbers are summarized above 
 gam=1.4  
 [rhoexact,uexact,pexact,eexact,xexact] = R.exact(gam,max_time[case],case,"rusanov"+str(case)) 

fig9, axs9 = plt.subplots(4,figsize=(6,8))
fig9.suptitle("Comparison of MacCormack FD and FV, Test 1, $Δx = 0.025$, $c_{max} = 1.0$",fontsize=10)

[rho1,u1,p1,e1,x1] = maccormackFD("Case 1",1.0,0.025,1.4)
[rho2,u2,p2,e2,x2] = maccormackFV("Case 1",1.0,0.025,1.4)

axs9[0].plot(xexact, rhoexact,"--k")
axs9[1].plot(xexact, uexact,"--k")
axs9[2].plot(xexact, pexact,"--k")
axs9[3].plot(xexact, eexact,"--k",label="Analytical")

axs9[0].plot(x1, rho1)
axs9[0].grid()
axs9[1].plot(x1, u1)
axs9[1].grid()
axs9[2].plot(x1, p1)
axs9[2].grid()
axs9[3].plot(x1, e1,label="MacCormack FD")
axs9[3].grid()

axs9[0].plot(x2, rho2,"-.")
axs9[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs9[0].set_xlim([0, 1])
axs9[1].plot(x2, u2,"-.")
axs9[1].set(ylabel="$u$ $(m/s)$")
axs9[1].set_xlim([0, 1])
axs9[1].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs9[2].plot(x2, p2,"-.")
axs9[2].set(ylabel="$p$ $(Pa)$")
axs9[2].set_xlim([0, 1])
axs9[2].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs9[3].plot(x2, e2,"-.",label="MacCormack FV")
axs9[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs9[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs9[3].set_xlim([0, 1])
axs9[3].yaxis.set_major_locator(ticker.MultipleLocator(0.2))

fig9.subplots_adjust(top=0.95)
fig9.savefig('HW4 Part 6_1.png', dpi=300)

##########

fig10, axs10 = plt.subplots(4,figsize=(6,8))
fig10.suptitle("Comparison of shock-capturing methods, Test 1, $Δx = 0.025$, $c_{max} = 1.0$",fontsize=10)

[rho1,u1,p1,e1,x1] = maccormackFD("Case 1",1.0,0.025,1.4)
[rho2,u2,p2,e2,x2] = lax_fredrich("Case 1",1.0,0.025,1.4)
[rho3,u3,p3,e3,x3] = rusanov("Case 1",1.0,0.025,1.4)
[rho4,u4,p4,e4,x4] = maccormackFCT("Case 1",1.0,0.025,1.4)

axs10[0].plot(xexact, rhoexact,"--k")
axs10[1].plot(xexact, uexact,"--k")
axs10[2].plot(xexact, pexact,"--k")
axs10[3].plot(xexact, eexact,"--k",label="Analytical")

axs10[0].plot(x1, rho1)
axs10[0].grid()
axs10[1].plot(x1, u1)
axs10[1].grid()
axs10[2].plot(x1, p1)
axs10[2].grid()
axs10[3].plot(x1, e1,label="MacCormack FD")
axs10[3].grid()

axs10[0].plot(x2, rho2)
axs10[1].plot(x2, u2)
axs10[2].plot(x2, p2)
axs10[3].plot(x2, e2,label="Lax-Fredrich")

axs10[0].plot(x3, rho3)
axs10[1].plot(x3, u3)
axs10[2].plot(x3, p3)
axs10[3].plot(x3, e3,label="Rusanov")

axs10[0].plot(x4, rho4)
axs10[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs10[0].set_xlim([0, 1])
axs10[1].plot(x4, u4)
axs10[1].set(ylabel="$u$ $(m/s)$")
axs10[1].set_xlim([0, 1])
axs10[1].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs10[2].plot(x4, p4)
axs10[2].set(ylabel="$p$ $(Pa)$")
axs10[2].set_xlim([0, 1])
axs10[2].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
axs10[3].plot(x4, e4,label="MacCormack FCT")
axs10[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs10[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=3)
axs10[3].set_xlim([0, 1])
axs10[3].yaxis.set_major_locator(ticker.MultipleLocator(0.2))

fig10.subplots_adjust(top=0.95)
fig10.savefig('HW4 Part 6_2.png', dpi=300)

##########