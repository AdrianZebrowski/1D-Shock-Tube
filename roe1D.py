def roe1D(test,c_max,dx):
    import numpy as np
    
    def roeavg(betaL,betaR,rhoL,rhoR):
        beta_tilde = (betaL*np.sqrt(rhoL)+betaR*np.sqrt(rhoR))/(np.sqrt(rhoL)+np.sqrt(rhoR))
        return beta_tilde
    
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
        
    while t <= t_final:
        for i in range(0,ix): # Loop for calculating Fn 
                
            im=i-1 
            ip=i+1
                
            # Neumann BCs implemented here using a series of if statements
            if i == 0:
                im = 1
            if i == ix-1:
                ip = ix-2
            
            rho_tilde = np.sqrt(rho[i]*rho[ip])
            u_tilde = roeavg(u[i],u[ip],rho[i],rho[ip])
            a_tilde = roeavg(a[i],a[ip],rho[i],rho[ip])
            H_tilde = roeavg(H[i],H[ip],rho[i],rho[ip])
            
            del_p = p[ip]-p[i]
            del_u = u[ip]-u[i]
            del_rho = rho[ip]-rho[i]
            
            alpha1_tilde = (del_p-a_tilde*rho_tilde*del_u)/(2*a_tilde**2)
            alpha2_tilde = (del_p+a_tilde*rho_tilde*del_u)/(2*a_tilde**2)
            alpha3_tilde = del_rho - del_p/a_tilde**2
            
            lambda1 = np.abs(u_tilde-a_tilde)
            lambda2 = np.abs(u_tilde+a_tilde)
            lambda3 = np.abs(u_tilde)
            
            eps1 = max(0,((u_tilde-a_tilde)-(u[i]-a[i])),((u[ip]-a[ip])-(u_tilde-a_tilde)))
            eps2 = max(0,((u_tilde+a_tilde)-(u[i]+a[i])),((u[ip]+a[ip])-(u_tilde+a_tilde)))
            
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