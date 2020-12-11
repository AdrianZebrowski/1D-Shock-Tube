def vanleer1DMUSCL(test,c_max,dx):
    import numpy as np

    if test == 1:
        # case 1 - Sod problem
        rhoL=1.0 
        uL=0.0 
        pL=1.0 
        rhoR=0.125 
        uR=0.0 
        pR=0.1 
        t_final = 0.25
    if test == 2:
        # case 2 - 123 problem - expansion left and expansion right
        rhoL=1.0 
        uL=-2.0
        pL=0.4
        rhoR=1.0
        uR=2.0
        pR=0.4
        t_final = 0.15
    if test == 3:
        # case 3 - blast problem - shock right, expansion left
        rhoL=1.0 
        uL=0.0
        pL=1000
        rhoR=1.0
        uR=0.
        pR=0.01
        t_final = 0.012
    if test == 4:
        # case 4 - blast problem - shock left, expansion right
        rhoL=1.0 
        uL=0.0
        pL=0.01
        rhoR=1.0
        uR=0.
        pR=100
        t_final = 0.035
    if test == 5:
        # case 5 - shock collision - shock left and shock right
        rhoL=5.99924
        uL=19.5975
        pL=460.894
        rhoR=5.99242
        uR=-6.19633
        pR=46.0950
        t_final = 0.035
        
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
    
    t = 0
    dt = c_max*dx/np.max(np.abs(u)+a)
    
    #calculate initial Q and F vectors
    Q[0,:] = rho
    Q[1,:] = rho*u
    Q[2,:] = rho*(e+u**2/2)
    
    rhoL = np.zeros(ix)
    rhoR = np.zeros(ix)
    uL = np.zeros(ix)
    uR = np.zeros(ix)
    pL = np.zeros(ix)
    pR = np.zeros(ix)
    
    rhoL[:] = rho
    rhoR[:] = rho
    uL[:] = u
    uR[:] = u
    pL[:] = p
    pR[:] = p
    
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
            
            rhoL[i],rhoR[i] = MUSCL(rho[im],rho[i],rho[ip],rho[ip2])
            uL[i],uR[i] = MUSCL(u[im],u[i],u[ip],u[ip2])  
            pL[i],pR[i] = MUSCL(p[im],p[i],p[ip],p[ip2])  
    
        aL = np.sqrt(gamma*pL/rhoL)
        aR = np.sqrt(gamma*pR/rhoR)
        ML = uL/aL
        MR = uR/aR
        
        for i in range(0,ix): # Loop for calculating Fn 
                    im=i-1 
                    ip=i+1
                        
                    # Neumann BCs implemented here using a series of if statements
                    if i == 0:
                        im = 1
                    if i == ix-1:
                        ip = ix-2
                                 
                    F_plus[0,i] = ((1/4)*rhoL[i]*aL[i]*(1+ML[i])**2)
                    F_plus[1,i] = ((1/4)*rhoL[i]*aL[i]*(1+ML[i])**2)*(2*aL[i]/gamma)*(((gamma-1)/2)*ML[i]+1)
                    F_plus[2,i] = ((1/4)*rhoL[i]*aL[i]*(1+ML[i])**2)*(2*aL[i]**2/(gamma**2-1))*(((gamma-1)/2)*ML[i]+1)**2
                    
                    F_minus[0,i] = -((1/4)*rhoR[i]*aR[i]*(1-MR[i])**2)
                    F_minus[1,i] = -((1/4)*rhoR[i]*aR[i]*(1-MR[i])**2)*(2*aR[i]/gamma)*(((gamma-1)/2)*MR[i]-1)
                    F_minus[2,i] = -((1/4)*rhoR[i]*aR[i]*(1-MR[i])**2)*(2*aR[i]**2/(gamma**2-1))*(((gamma-1)/2)*MR[i]-1)**2
                
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
        dt = c_max*dx/np.max(np.abs(u)+a)
        t = t+dt
        
        print(t)
        
    return rho, u, p, e, x