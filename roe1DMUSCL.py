def roe1DMUSCL(test,c_max,dx):
    import numpy as np
    
    def roeavg(betaL,betaR,rhoL,rhoR):
        beta_tilde = (betaL*np.sqrt(rhoL)+betaR*np.sqrt(rhoR))/(np.sqrt(rhoL)+np.sqrt(rhoR))
        return beta_tilde
    
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
    F = np.zeros((3,ix))
    F_half = np.zeros((3,ix))
    
    #establish initial condition vectors
    rho[0:half] = rhoL
    rho[half:ix] = rhoR
    u[0:half] = uL
    u[half:ix] = uR
    p[0:half] = pL
    p[half:ix] = pR
    e = p/((gamma-1)*rho)
    H = e+p/rho+u**2/2
    a = np.sqrt((gamma*p/rho))
    
    t = 0
    dt = c_max*dx/np.max(np.abs(u)+a)
    
    #calculate initial Q and F vectors
    Q[0,:] = rho
    Q[1,:] = rho*u
    Q[2,:] = rho*(e+u**2/2)
    
    F[0,:] = rho*u
    F[1,:] = rho*u**2+p
    F[2,:] = rho*u*(e+p/rho+u**2/2)
    
    rhoL = np.zeros(ix)
    rhoR = np.zeros(ix)
    uL = np.zeros(ix)
    uR = np.zeros(ix)
    pL = np.zeros(ix)
    pR = np.zeros(ix)
    eL = np.zeros(ix)
    eR = np.zeros(ix)
    
    rhoL[:] = rho
    rhoR[:] = rho
    uL[:] = u
    uR[:] = u
    pL[:] = p
    pR[:] = p
    eL[:] = e
    eR[:] = e
        
    while t <= t_final:
        for i in range(1,ix-2): # Loop for calculating Fn  
            
            im=i-1 
            ip=i+1
            ip2=i+2
            
            rhoL[i],rhoR[i] = MUSCL(rho[im],rho[i],rho[ip],rho[ip2])
            uL[i],uR[i] = MUSCL(u[im],u[i],u[ip],u[ip2])  
            pL[i],pR[i] = MUSCL(p[im],p[i],p[ip],p[ip2])  
            eL[i],eR[i] = MUSCL(e[im],e[i],e[ip],e[ip2])
        
        aL = np.sqrt(gamma*pL/rhoL)
        aR = np.sqrt(gamma*pR/rhoR)
        HL = eL+pL/rhoL+uL**2/2
        HR = eR+pR/rhoR+uR**2/2
        
        for i in range(0,ix): # Loop for calculating Fn 
                
            im=i-1 
            ip=i+1
                
            # Neumann BCs implemented here using a series of if statements
            if i == 0:
                im = 1
            if i == ix-1:
                ip = ix-2
            
            rho_tilde = np.sqrt(rhoL[i]*rhoR[i])
            u_tilde = roeavg(uL[i],uR[i],rhoL[i],rhoR[i])
            a_tilde = roeavg(aL[i],aR[i],rhoL[i],rhoR[i])
            H_tilde = roeavg(HL[i],HR[i],rhoL[i],rhoR[i])
            
            del_p = pR[i]-pL[i]
            del_u = uR[i]-uL[i]
            del_rho = rhoR[i]-rhoL[i]
            
            alpha1_tilde = (del_p-a_tilde*rho_tilde*del_u)/(2*a_tilde**2)
            alpha2_tilde = (del_p+a_tilde*rho_tilde*del_u)/(2*a_tilde**2)
            alpha3_tilde = del_rho - del_p/a_tilde**2
            
            lambda1 = np.abs(u_tilde-a_tilde)
            lambda2 = np.abs(u_tilde+a_tilde)
            lambda3 = np.abs(u_tilde)
            
            eps1 = max(0,((u_tilde-a_tilde)-(uL[i]-aL[i])),((uR[i]-aR[i])-(u_tilde-a_tilde)))
            eps2 = max(0,((u_tilde+a_tilde)-(uL[i]+aL[i])),((uR[i]+aR[i])-(u_tilde+a_tilde)))
            
            if eps1 > lambda1:
                lambda1 = eps1
            if eps2 > lambda2:
                lambda2 = eps2
            
            F_half[0,i] = 0.5*(F[0,i]+F[0,ip]) - 0.5*(alpha1_tilde*lambda1+alpha2_tilde*lambda2+alpha3_tilde*lambda3)
            F_half[1,i] = 0.5*(F[1,i]+F[1,ip]) - 0.5*(alpha1_tilde*lambda1*(u_tilde-a_tilde)+alpha2_tilde*lambda2*(u_tilde+a_tilde)+alpha3_tilde*lambda3*u_tilde)
            F_half[2,i] = 0.5*(F[2,i]+F[2,ip]) - 0.5*(alpha1_tilde*lambda1*(H_tilde-u_tilde*a_tilde)+alpha2_tilde*lambda2*(H_tilde+u_tilde*a_tilde)+alpha3_tilde*lambda3*(u_tilde**2/2))
            
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
        H = e+p/rho+u**2/2
        a = np.sqrt((gamma*p/rho))
        
        F[0,:] = rho*u
        F[1,:] = rho*u**2+p
        F[2,:] = rho*u*(e+p/rho+u**2/2)
        
        dt = c_max*dx/np.max(np.abs(u)+a)
        t = t+dt
        
        print(t)
        
    return rho, u, p, e, x