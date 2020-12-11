def vanleer1D(test,c_max,dx):
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
    M = u/a
    
    t = 0
    dt = c_max*dx/np.max(np.abs(u)+a)
    
    #calculate initial Q and F vectors
    Q[0,:] = rho
    Q[1,:] = rho*u
    Q[2,:] = rho*(e+u**2/2)
        
    while t <= t_final:
        for i in range(0,ix): # Loop for calculating Fn 
                
            im=i-1 
            ip=i+1
                
            # Neumann BCs implemented here using a series of if statements
            if i == 0:
                im = 1
            if i == ix-1:
                ip = ix-2
                    
            F_plus[0,i] = ((1/4)*rho[i]*a[i]*(1+M[i])**2)
            F_plus[1,i] = ((1/4)*rho[i]*a[i]*(1+M[i])**2)*(2*a[i]/gamma)*(((gamma-1)/2)*M[i]+1)
            F_plus[2,i] = ((1/4)*rho[i]*a[i]*(1+M[i])**2)*(2*a[i]**2/(gamma**2-1))*(((gamma-1)/2)*M[i]+1)**2
            
            F_minus[0,i] = -((1/4)*rho[ip]*a[ip]*(1-M[ip])**2)
            F_minus[1,i] = -((1/4)*rho[ip]*a[ip]*(1-M[ip])**2)*(2*a[ip]/gamma)*(((gamma-1)/2)*M[ip]-1)
            F_minus[2,i] = -((1/4)*rho[ip]*a[ip]*(1-M[ip])**2)*(2*a[ip]**2/(gamma**2-1))*(((gamma-1)/2)*M[ip]-1)**2
        
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
        
    return rho, u, p, e, x