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
 
rho1, u1, p1, e1, x1 = vanleer1DMUSCL(1,0.3,0.01)  
rho2, u2, p2, e2, x2 = roe1DMUSCL(1,0.3,0.01)

fig1, axs1 = plt.subplots(4,figsize=(6,8))
fig1.suptitle("Comparison of shock-capturing methods, Test 1, $Δx = 0.025$, $c_{max} = 0.5$",fontsize=10)

axs1[0].plot(xexact, rhoexact,"--k")
axs1[1].plot(xexact, uexact,"--k")
axs1[2].plot(xexact, pexact,"--k")
axs1[3].plot(xexact, eexact,"--k",label="Analytical")

axs1[0].plot(x1, rho1)
axs1[1].plot(x1, u1)
axs1[2].plot(x1, p1)
axs1[3].plot(x1, e1,label="Van Leer")

axs1[0].plot(x1, rho2,'-.')
axs1[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs1[0].set_xlim([0, 1])
axs1[0].grid()
axs1[1].plot(x1, u2,'-.')
axs1[1].set(ylabel="$u$ $(m/s)$")
axs1[1].set_xlim([0, 1])
axs1[1].grid()
axs1[2].plot(x1, p2,'-.')
axs1[2].set(ylabel="$p$ $(Pa)$")
axs1[2].set_xlim([0, 1])
axs1[2].grid()
axs1[3].plot(x1, e2,'-.',label="Roe")
axs1[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs1[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs1[3].set_xlim([0, 1])
axs1[3].grid()

fig1.subplots_adjust(top=0.95)    

rho1, u1, p1, e1, x1 = vanleer1D(5,0.5,0.0125)  
rho2, u2, p2, e2, x2 = roe1D(5,0.5,0.0125)

fig2, axs1 = plt.subplots(4,figsize=(6,8))
fig2.suptitle("Comparison of shock-capturing methods, Test 1, $Δx = 0.0125$, $c_{max} = 0.5$",fontsize=10)

axs1[0].plot(xexact, rhoexact,"--k")
axs1[1].plot(xexact, uexact,"--k")
axs1[2].plot(xexact, pexact,"--k")
axs1[3].plot(xexact, eexact,"--k",label="Analytical")

axs1[0].plot(x1, rho1)
axs1[1].plot(x1, u1)
axs1[2].plot(x1, p1)
axs1[3].plot(x1, e1,label="Van Leer")

axs1[0].plot(x1, rho2,'-.')
axs1[0].set(ylabel="$ρ$ $(kg/m^3)$")
axs1[0].set_xlim([0, 1])
axs1[0].grid()
axs1[1].plot(x1, u2,'-.')
axs1[1].set(ylabel="$u$ $(m/s)$")
axs1[1].set_xlim([0, 1])
axs1[1].grid()
axs1[2].plot(x1, p2,'-.')
axs1[2].set(ylabel="$p$ $(Pa)$")
axs1[2].set_xlim([0, 1])
axs1[2].grid()
axs1[3].plot(x1, e2,'-.',label="Roe")
axs1[3].set(xlabel="$x$", ylabel="$e$ $(J/kg∙K)$")
axs1[3].legend(loc='upper center', bbox_to_anchor=(0.5,-0.3),ncol=4)
axs1[3].set_xlim([0, 1])
axs1[3].grid()

fig2.subplots_adjust(top=0.95)    