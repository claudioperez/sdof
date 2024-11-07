"""
Copyright (c) 2018
    George Papazafeiropoulos
    Captain, Infrastructure Engineer, Hellenic Air Force
    Civil Engineer, M.Sc., Ph.D. candidate, NTUA
    Email: gpapazafeiropoulos@yahoo.gr

https://www.mathworks.com/matlabcentral/fileexchange/44649-general-single-step-single-solve-integration-algorithm
"""
import numpy as np
import warnings

class GSSSS:
    def __init__(self,
                 parameters,
                 predictor=None, 
                 corrector=None,
                 scheme=None,
                 maxtol=1e-6,
                 maxiter=10):
        
        if isinstance(parameters, float):
            rinf = parameters 

        elif len(parameters) == 2:
            beta, gamma = parameters
        


def gssss(dt: float,
          force,
          model,
          m,k,
          u0=0.0, 
          v0=0.0,
          maxtol=1e-6,
          maxiter=10):
    """
    Non Linear Implicit Dynamic Analysis of a bilinear kinematic hardening
    hysteretic structure with elastic damping

       [u,ut,utt,Fs,Ey,Es,Ed,jiter] = NLIDABLKIN(dt,xgtt,m,k_hi,k_lo,uy,ksi,...
           AlgID,u0,ut0,maxtol,maxiter,dak)
           General linear implicit direct time integration of second order
           differential equations of a bilinear elastoplastic hysteretic SDOF
           dynamic system with elastic damping, with lumped mass.

       Description
           The General Single Step Single Solve (GSSSS) family of algorithms
           published by X.Zhou & K.K.Tamma (2004) is employed for direct time
           integration of the general linear or nonlinear structural Single
           Degree of Freedom (SDOF) dynamic problem. Selection among 9
           algorithms, all designed according to the above journal article, can
           be made in this routine. These algorithms encompass the scope of
           Linear Multi-Step (LMS) methods and are limited by the Dahlquist
           barrier theorem (Dahlquist, 1963).

       Input parameters
           dt (scalar) is the time step of the integration
           xgtt ([#NumSteps# x 1]) is the acceleration time history which is
               imposed at the lumped mass of the SDOF structure.
           #AlgID# (string as follows) is the algorithm to be used for the time
               integration. It can be one of the following strings for superior
               optimally designed algorithms:
                   'generalized a-method': The generalized a-method (Chung & Hulbert, 1993)
                   'HHT a-method': The Hilber-Hughes-Taylor method (Hilber, Hughes & Taylor, 1977)
                   'WBZ': The Wood-Bossak-Zienkiewicz method (Wood, Bossak & Zienkiewicz, 1980)
                   'U0-V0-Opt': Optimal numerical dissipation and dispersion zero order displacement zero order velocity algorithm
                   'U0-V0-CA': Continuous acceleration (zero spurious root at
                   the low frequency limit) zero order displacement zero order
                   velocity algorithm
                   'U0-V0-DA': Discontinuous acceleration (zero spurious root at
                               the high frequency limit) zero order displacement zero order
                               velocity algorithm
                   'U0-V1-Opt': Optimal numerical dissipation and dispersion
                                zero order displacement first order velocity algorithm
                   'U0-V1-CA': Continuous acceleration (zero spurious root at
                                the low frequency limit) zero order displacement first order
                                velocity algorithm
                   'U0-V1-DA': Discontinuous acceleration (zero spurious root at
                                the high frequency limit) zero order displacement first order
                                velocity algorithm
                   'U1-V0-Opt': Optimal numerical dissipation and dispersion
                                first order displacement zero order velocity algorithm
                   'U1-V0-CA': Continuous acceleration (zero spurious root at
                   the low frequency limit) first order displacement zero order
                   velocity algorithm
                   'U1-V0-DA': Discontinuous acceleration (zero spurious root at
                   the high frequency limit) first order displacement zero order
                   velocity algorithm
                   'Newmark ACA': Newmark Average Constant Acceleration method
                   'Newmark LA': Newmark Linear Acceleration method
                   'Newmark BA': Newmark Backward Acceleration method
                   'Fox-Goodwin': Fox-Goodwin formula
           #u0# (scalar) is the initial displacement.
           #v0# (scalar) is the initial velocity.
           #rinf# (scalar) is the minimum absolute value of the eigenvalues of
               the amplification matrix. For the amplification matrix see
               eq.(61) in Zhou & Tamma (2004).
           #maxtol# (scalar) is the maximum tolerance of convergence of the Full
               Newton Raphson method for numerical computation of acceleration.
           #maxiter# (scalar) is the maximum number of iterations per increment. If
               #maxiter#=0 then iterations are not performed and the #maxtol#
               parameter is not taken into account.
           #dak# (scalar) is the infinitesimal acceleration for the
               calculation of the derivetive required for the convergence of the
               Newton-Raphson iteration.

       Output parameters
           #u# ([1 x #NumSteps#]) is the time-history of displacement
           #ut# ([1 x #NumSteps#]) is the time-history of velocity
           #utt# ([1 x #NumSteps#]) is the time-history of acceleration
           #Fs# ([1 x #NumSteps#]) is the time-history of the internal
               force of the structure analysed.
           #Ey# ([1 x #NumSteps#]) is the time history of the sum of the
               energy dissipated by yielding during each time step and the
               recoverable strain energy of the system (incremental).
               cumsum(#Ey#)-#Es# gives the time history of the total energy
               dissipated by yielding from the start of the dynamic analysis.
           #Es# ([1 x #NumSteps#]) is the time-history of the recoverable
               strain energy of the system (total and not incremental).
           #Ed# ([1 x #NumSteps#]) is the time-history of the energy
               dissipated by viscoelastic damping during each time step
               (incremental). cumsum(#Ed#) gives the time history of the total
               energy dissipated from the start of the dynamic analysis.
           #jiter# ([1 x #NumSteps#]) is the iterations per increment

       Notation in the code
           u   =displacement
           un  =displacement after increment n
           ut  =velocity
           vn =velocity after increment n
           utt =acceleration
           an=acceleration after increment n

      __________________________________________________________________________
    _________________________________________________________________________"""
    if isinstance(k, tuple):
        k_hi,k_lo,ksi = k

    # required inputs
        # if ~isscalar(dt)
        #     error('dt is not scalar')
        #
        # if dt<=0
        #     error('dt is zero or negative')
        #
        # if ~isvector(xgtt)
        #     error('xgtt is not vector')
        #
        # if ~isscalar(m)
        #     error('m is not scalar')
        #
        # if m<=0
        #     error('m is zero or negative')
        #
        # if rinf<0 || rinf>1
        #     error('rinf is lower than 0 or higher than 1')
        #
        # if ~isscalar(maxtol)
        #     error('maxtol is not scalar')
        #
        # if maxtol<=0
        #     error('maxtol is zero or negative')
        #
        # if ~isscalar(dak)
        #     error('dak is not scalar')
        #
        # if dak<=0
        #     error('dak is zero or negative')
        #

    alg  = ""
    rinf = 1.0
    W1, W1L1, W2L2, W3L3, W1L4, W2L5, W1L6, l1, l2, l3, l4, l5 = parameters(alg, rinf)


    ## Calculation
    # number of analysis increments
    nt = len(force)
    # Initialize output
    u     = np.zeros((1,nt))
    v     = np.zeros((1,nt))
    a     = np.zeros((1,nt))
    Fs    = np.zeros((1,nt))
    Ey    = np.zeros((1,nt))
    Es    = np.zeros((1,nt))
    Ed    = np.zeros((1,nt))
    # set initial values of displacement, velocity, acceleration (u0,ut0 and
    # utt0 respectively) at n=0.
    u[0]  = u0
    v[0]  = v0

    # construct lumped mass matrix
    M = m; # diag(m,0)
    state = 0
    pn,Kn,Cn,k,state = model(u0,v0,state)
    a[0] = -force[0] - np.linalg.solve(M, pn)

    Fs[0] = pn
    # initial assignments
    uo   = u0
    vo   = v0
    ao   = a[0]

    # Begin integration
    for n in range(1, nt):
        # effective force
        g = -pn \
                - Kn*(W1L1*dt*vo + W2L2*dt**2*ao) \
                - W1L4*dt*Cn*ao \
                + M*(((1-W1)*force[n] + W1*force[n+1]) - ao)
        # effective mass
        Dg = W1L6*M + W2L5*dt*Cn + W3L3*dt**2*Kn
        # initial estimate of da
        dan = np.linalg.solve(Dg, g)

        # start iteration number
        j=1
        # set initial quotient of variation of da equal to maxtol
        quda = maxtol

        # full Newton-Raphson iterations
        while abs(quda) >= maxtol and j<=maxiter:
            # _________________________________________________________________
            #/
            # State and residual at an+dan
            un = uo + vo*l1*dt + l2*ao*dt**2 + l3*dan*dt**2
            vn = vo + ao*l4*dt + l5*dan*dt
            an = ao + dan
            # force due to stiffness and damping
            pn1k, Kn1k, Cn1k,_,_ = model(un,vn,state)
            # effective force
            g = -pn1k \
                -Kn1k*(W1L1*dt*vn + W2L2*dt**2*an) \
                -Cn1k*W1L4*dt*an \
                +M*((1-W1)*force[n] + W1*force[n+1]-an)
            # effective mass
            Meffn1k = Kn1k*W3L3*dt**2 + Cn1k*W2L5*dt + M*W1L6
            # residual
            Rn1k = g - Meffn1k*dan
            #\_________________________________________________________________

            # _________________________________________________________________
            #/
            # displacement, velocity, acceleration, internal force, stiffness
            # and damping for #an+(dan+dak)#
            # calculate the derivative at #an+dan# as:
            # #dR/da=(dRn1k-Rn1k)/(an+(dan+dak)-(an+dan))=(dRn1k-Rn1k)/dak#
            # update kinematic quantities
            dun1k = uo + l1*vo*dt + l2*ao*dt**2 + l3*(dan + dak)*dt**2
            dvn1k = vo + l4*ao*dt + l5*(dan + dak)*dt
            dan1k = ao + (dan + dak)

            # force due to stiffness and damping
            dpn1k,dKn1k,dCn1k,_, __ = model(dun1k,dvn1k,state)
            # effective force
            dFeffn1k = -dpn1k \
                - dKn1k*(W1L1*dt*dvn1k + W2L2*dt**2*dan1k) \
                - dCn1k*W1L4*dt*dan1k \
                + M*((1-W1)*force[n] + W1*force[n+1] - dan1k)
            # effective mass
            dMeffn1k = dKn1k*W3L3*dt**2 + dCn1k*W2L5*dt + M*W1L6
            # residual
            dRn1k = dFeffn1k - dMeffn1k*dan1k
            #\_________________________________________________________________

            # Full Newton-Raphson update:
            # #da_new = da - Rn1k/(dR/da)=da*(1-Rn1k/(dRn1k/dak)/da)#
            # (to be checked for while loop termination)
            quda = (Rn1k/(dRn1k - Rn1k)*dak)/dan

            # test if derivative becomes zero
            if any(np.isinf(quda)):
                break
                #quda = np.zeros(size(quda))

            # update da
            dan = (1 - quda) * dan
            # update iteration number
            j += 1

        # _____________________________________________________________________
        #/
        # displacement and its derivatives after iteration k+1 of increment
        # n+1
        un = uo + l1*vo*dt + l2*ao*dt**2 + l3*dan*dt**2
        vn = vo + l4*ao*dt + l5*dan*dt
        an = ao + dan
        # internal force, stiffness and damping after iteration k+1 of
        # increment n+1
        pn1k,Kn1k,Cn1k,k,state = model(un,vn,state)
        # _____________________________________________________________________

        # assignments to output parameters
        u[n+1] = un
        v[n+1] = vn
        a[n+1] = an
        Fs[n+1]=pn1k
        if False:
            Ey[n+1] =-(np.cumsum(pn1k-Cn1k*vn) + np.cumsum(pn-Cn*vo))/2.0*diff([[un-uo],[0]])
            Es[n+1] =  np.cumsum(pn1k-Cn1k*vn)**2 / k_hi/2
            Ed[n+1] =-(np.cumsum(Cn1k*vn)+np.cumsum(Cn*vo))/2.0*diff([[un-uo],[0]])

        # assignments for next increment
        pn=pn1k
        Kn=Kn1k
        Cn=Cn1k
        uo=un
        vo=vn
        ao=an


def model(u,ut,k_hi,k_lo,uy,M,ksi,k,d):
    """Bilinear elastoplastic hysteretic model with elastic viscous damping

    [f,K,C,k,d] = BLKIN(u,ut,k_hi,k_lo,uy,M,ksi,k,d)
        Define the internal force vector, tangent stiffness matrix and
        tangent damping matrix of a bilinear elastoplastic hysteretic
        structure with elastic damping as a function of displacement and
        velocity.

    Description
    ---------------------------
        The MDOF structure modeled with this function consists of lumped
        masses connected with stiffness and damping elements in series. Each
        lumped mass has one degree of freedom. The first degree of freedom is
        at the top of the structure and the last at its fixed base. However,
        the last degree of freedom is not included in the input arguments of
        the function, i.e. not contained in #ndof#, as it is always fixed.
        The nonlinear stiffness is virtually of the bilinear type, where an
        initial stiffness and a post-yield stiffness are defined. The
        unloading or reloading curve of this model are parallel to the
        initial loading curve, and a hysteresis loop is created by
        continuously loading and unloading the structure above its yield
        limit. This behavior can be viewed as hardening of the kinematic
        type.
        An appropriate reference for this function definition is Hughes,
        Pister & Taylor (1979): "Implicit-explicit finite elements in
        nonlinear transient analysis". This function should be defined in
        accordance with equations (3.1), (3.2) and (3.3) of this paper. This
        representation has as special cases nonlinear elasticity and a class
        of nonlinear rate-type viscoelastic materials. Tangent stiffness and
        tangent damping matrices are the "consistent" linearized operators
        associated to f in the sense of [Hughes & Pister, "Consistent
        linearization in mechanics of solids", Computers and Structures, 8
        (1978) 391-397].

           m (scalar) is the lumped masses of the structure. Define the
               lumped masses from the top to the bottom, excluding the fixed dof
               at the base
           #k_hi# (scalar): is the initial stiffness of the system before
               its first yield, i.e. the high stiffness. Give the stiffness of
               each storey from top to bottom.
           #k_lo# (scalar): is the post-yield stiffness of the system,
               i.e. the low stiffness. Give the stiffness of each storey from
               top to bottom.
           #uy# (scalar): is the yield limit of the stiffness elements of
               the structure. The element is considered to yield, if the
               interstorey drift between degrees of freedom i and i+1 exceeds
               uy(i). Give the yield limit of each storey from top to bottom.
           #ksi# (scalar): ratio of critical viscous damping of the system,
               assumed to be unique for all damping elements of the structure.

    Input parameters
    ---------------------
        #u# (scalar): absolute displacement.
        #ut# (scalar): absolute velocity.
        #k_hi# (scalar): initial stiffness of the system before its first
            yield, i.e. the high stiffness.
        #k_lo# (scalar): post-yield stiffness of the system, i.e. the low
            stiffness.
        #uy# (scalar): yield limit of the structure. The structure is
            considered to yield, if the displacement exceeds uy(i).
        #M# (scalar): lumped mass.
        #ksi# (scalar): ratio of critical viscous damping of the system,
            assumed to be unique for all damping elements of the structure.
        #k# (scalar): is the stiffness vector which takes into account
            any plastic response of the structure. It is used to record the
            status of the structure so that it is known before the next
            application of this function at a next (time) step. Initialize by
            setting #k#=#k_hi#.
        #d# (scalar): is the equilibrium displacement vector which takes into
            account any plastic response of the structure. It is used to
            record the status of the structure so that it is known before the
            next application of this function at a next (time) step.
            Initialize by setting #d#=zeros(#ndof#,1).

    Output parameters
    ----------
        #f# (scalar): internal force vector of the structure (sum of forces
            due to stiffness and damping) at displacement #u# and velocity
            #ut#
        #K# (scalar): tangent stiffness matrix (nonlinear function of
            displacement #u# and velocity #ut#). It is equivalent to the
            derivative d(#f#)/d(#u#)
        #C# (scalar): tangent damping matrix (nonlinear function of
            displacement #u# and velocity #ut#). It is equivalent to the
            derivative d(#f#)/d(#u#)
        #k# (scalar): is the stiffness vector which takes into account any
            plastic response of the structure. It is used to record the
            status of the structure so that it is known before the next
            application of this function at a next (time) step.
        #d# (scalar): is the equilibrium displacement vector which takes into
            account any plastic response of the structure. It is used to
            record the status of the structure so that it is known before the
            next application of this function at a next (time) step.

    Verification:
    -------------------------
        u=0:0.2:4
        ut=0.001*ones(1,numel(u))
        u=[u,u(end:-1:1)]
        ut=[ut,-ut]
        u=[u,-u]
        ut=[ut,ut(end:-1:1)]
        u=[u u]
        ut=[ut ut]
        k_hi=1000
        k_lo=1
        uy=2
        M=1
        ksi=0.05
        k=k_hi
        d=0
        f=zeros(1,numel(u))
        for i=1:numel(u)
            [f(i),K,C,k,d] = BLKIN(u(i),ut(i),k_hi,k_lo,uy,M,ksi,k,d)

        figure()
        plot(u,f)

    _________________________________________________________________________"""
        # if ~isscalar(k_hi)
        #     error('k_hi is not scalar')
        #
        # if k_hi<=0
        #     error('k_hi is zero or negative')
        #
        # if ~isscalar(k_lo)
        #     error('k_lo is not scalar')
        #
        # if k_lo<=0
        #     error('k_lo is zero or negative')
        #
        # if ~isscalar(uy)
        #     error('uy is not scalar')
        #
        # if uy<=0
        #     error('uy is zero or negative')
        #
        # # optional inputs
        # if ~isscalar(ksi)
        #     error('ksi is not scalar')
        #
        # if ksi<=0
        #     error('ksi is zero or negative')
        #

    # Elastic tangent stiffness matrix
    K=k_hi
    # Elastic tangent damping matrix
    C=2*ksi* np.sqrt(K*M)
    # force from stiffness (not damping) of the current storey
    fK=k*(u(1)-d)
    # eq.(46) in ...
    fy=k_lo*(u(1))+(k_hi-k_lo)*(uy*np.sign(ut(1)))
    # check for yielding or load reversal
    if k == k_hi and ut(1)>0 and fK>fy:
        # the system has just exceeded its positive yield force level
        k = k_lo
        d = (1-k_hi/k_lo)*uy
    elif k == k_hi and ut(1)<0 and fK<fy:
        # check for yielding
        # the system has just exceeded its negative yield force level
        k=k_lo
        d=(k_hi/k_lo-1)*uy
    elif k==k_lo and fK*(ut(1))<0:
        # check for load reversal
        # the system reloads from negative ultimate displacement or unloads
        # from positive ultimate displacement
        k=k_hi
        d=(u(1))-k_lo/k_hi*(u(1)-d)

    fK_bak=k*(u(1)-d)
    # Update the elastic tangent stiffness matrix
    K=k
    # internal force due to stiffness and damping
    f = fK_bak + C*ut
    return f,K,C,k,d


def parameters(alg, rinf):
    """
    """
    if alg == "U0-V0-Opt":
        # zero-order displacement & velocity overshooting behavior and
        # optimal numerical dissipation and dispersion
        # rinf must belong to [0 1]
        if rinf<0:
            rinf=0
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 0')

        if rinf>1:
            rinf=1; # mid-point rule a-form algorithm
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1')

        w1 = -15*(1-2*rinf)/(1-4*rinf); # suggested
        w2 = 15*(3-4*rinf)/(1-4*rinf); # suggested
        w3 = -35*(1-rinf)/(1-4*rinf); # suggested
        W1L1 = 1/(1+rinf)
        W2L2 = 1/2/(1+rinf)
        W3L3 = 1/2/(1+rinf)**2
        W1L4 = 1/(1+rinf)
        W2L5 = 1/(1+rinf)**2; # suggested
        W1L6=(3-rinf)/2/(1+rinf)
        l1 = 1
        l2 = 1/2
        l3 = 1/2/(1+rinf)
        l4 = 1
        l5 = 1/(1+rinf)

    elif alg == "U0-V0-CA":
        # zero-order displacement & velocity overshooting behavior and
        # continuous acceleration
        # rinf must belong to [1/3 1]
        if rinf<1/3:
            rinf=1/3
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1/3')

        if rinf>1:
            rinf=1; # Newmark average acceleration a-form algorithm
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1')

        w1=-15*(1-5*rinf)/(3-7*rinf); # suggested
        w2 = 15*(1-13*rinf)/(3-7*rinf); # suggested
        w3 = 140*rinf/(3-7*rinf); # suggested
        W1L1=(1+3*rinf)/2/(1+rinf)
        W2L2=(1+3*rinf)/4/(1+rinf)
        W3L3=(1+3*rinf)/4/(1+rinf)**2
        W1L4=(1+3*rinf)/2/(1+rinf)
        W2L5=(1+3*rinf)/2/(1+rinf)**2; # suggested
        W1L6 = 1
        l1 = 1
        l2 = 1/2
        l3 = 1/2/(1+rinf)
        l4 = 1
        l5 = 1/(1+rinf)

    elif alg == "U0-V0-DA":
        # zero-order displacement & velocity overshooting behavior and
        # discontinuous acceleration
        # rinf must belong to [0 1]
        if rinf<0:
            rinf=0
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 0')

        if rinf>1:
            rinf=1; # Newmark average acceleration a-form algorithm
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1')

        w1=-15; # suggested
        w2 = 45; # suggested
        w3=-35; # suggested
        W1L1 = 1
        W2L2 = 1/2
        W3L3 = 1/2/(1+rinf)
        W1L4 = 1
        W2L5 = 1/(1+rinf); # suggested
        W1L6=(3+rinf)/2/(1+rinf)
        l1 = 1
        l2 = 1/2
        l3 = 1/2/(1+rinf)
        l4 = 1
        l5 = 1/(1+rinf)


    elif alg == "U0-V1-Opt":
        # zero-order displacement & first-order velocity overshooting
        # behavior and optimal numerical dissipation and dispersion
        # This is the generalized a-method (Chung & Hulbert, 1993)
        # rinf must belong to [0 1]
        if rinf<0:
            rinf=0
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 0')

        if rinf>1:
            rinf=1
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1')

        w1 = -15*(1-2*rinf)/(1-4*rinf)
        w2 =  15*(3-4*rinf)/(1-4*rinf)
        w3 = -35*(1-rinf)/(1-4*rinf)

        W1L1 = 1/(1+rinf)
        W2L2 = 1/2/(1+rinf)
        W3L3 = 1/(1+rinf)**3
        W1L4 = 1/(1+rinf)
        W2L5=(3-rinf)/2/(1+rinf)**2
        W1L6=(2-rinf)/(1+rinf)
        l1 = 1
        l2 = 1/2
        l3 = 1/(1+rinf)**2
        l4 = 1
        l5=(3-rinf)/2/(1+rinf)


    elif alg == "generalized a-method":
        # zero-order displacement & first-order velocity overshooting
        # behavior and optimal numerical dissipation and dispersion
        # This is the generalized a-method (Chung & Hulbert, 1993)
        # rinf must belong to [0 1]
        if rinf<0:
            rinf=0
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 0')

        if rinf>1:
            rinf=1
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1')

        w1 = -15*(1-2*rinf)/(1-4*rinf)
        w2 = 15*(3-4*rinf)/(1-4*rinf)
        w3 = -35*(1-rinf)/(1-4*rinf)
        W1L1 = 1/(1+rinf)
        W2L2 = 1/2/(1+rinf)
        W3L3 = 1/(1+rinf)**3
        W1L4 = 1/(1+rinf)
        W2L5=(3-rinf)/2/(1+rinf)**2
        W1L6=(2-rinf)/(1+rinf)
        l1 = 1
        l2 = 1/2
        l3 = 1/(1+rinf)**2
        l4 = 1
        l5=(3-rinf)/2/(1+rinf)


    elif alg == "U0-V1-CA":
        # zero-order displacement & first-order velocity overshooting
        # behavior and continuous acceleration
        # This is the Hilber-Hughes-Taylor method (Hilber, Hughes &
        # Taylor, 1977)
        # rinf must belong to [1/2 1]
        if rinf<1/2:
            rinf=1/2
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1/2')

        if rinf>1:
            rinf=1
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1')

        w1=-15*(1-2*rinf)/(2-3*rinf)
        w2 = 15*(2-5*rinf)/(2-3*rinf)
        w3=-35*(1-3*rinf)/2/(2-3*rinf)
        W1L1 = 2*rinf/(1+rinf)
        W2L2=rinf/(1+rinf)
        W3L3 = 2*rinf/(1+rinf)**3
        W1L4 = 2*rinf/(1+rinf)
        W2L5=rinf*(3-rinf)/(1+rinf)**2
        W1L6 = 1
        l1 = 1
        l2 = 1/2
        l3 = 1/(1+rinf)**2
        l4 = 1
        l5=(3-rinf)/2/(1+rinf)

    elif alg == "HHT a-method":
        # zero-order displacement & first-order velocity overshooting
        # behavior and continuous acceleration
        # This is the Hilber-Hughes-Taylor method (Hilber, Hughes &
        # Taylor, 1977)
        # rinf must belong to [1/2 1]
        if rinf<1/2:
            rinf=1/2
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1/2')

        if rinf>1:
            rinf=1
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1')

        w1   = -15*(1-2*rinf)/(2-3*rinf)
        w2   =  15*(2-5*rinf)/(2-3*rinf)
        w3   = -35*(1-3*rinf)/2/(2-3*rinf)
        W1L1 = 2*rinf/(1+rinf)
        W2L2 = rinf/(1+rinf)
        W3L3 = 2*rinf/(1+rinf)**3
        W1L4 = 2*rinf/(1+rinf)
        W2L5 = rinf*(3-rinf)/(1+rinf)**2
        W1L6 = 1.0
        l1   = 1.0
        l2   = 1/2
        l3   = 1/(1+rinf)**2
        l4   = 1.0
        l5   = (3-rinf)/2/(1+rinf)


    elif alg == "U0-V1-DA":
        # zero-order displacement & first-order velocity overshooting
        # behavior and discontinuous acceleration
        # This is the Wood烹ossak忙ienkiewicz method (Wood, Bossak &
        # Zienkiewicz, 1980)
        # rinf must belong to [0 1]
        if rinf<0:
            rinf=0
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 0')

        if rinf>1:
            rinf=1
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1')

        w1   = -15
        w2   =  45
        w3   = -35
        W1L1 = 1.0
        W2L2 = 1/2
        W3L3 = 1/(1+rinf)**2
        W1L4 = 1
        W2L5 = (3-rinf)/2/(1+rinf)
        W1L6 = 2/(1+rinf)
        l1   = 1.0
        l2   = 1/2
        l3   = 1/(1+rinf)**2
        l4   = 1.0
        l5   = (3-rinf)/2/(1+rinf)

    elif alg == "WBZ":
        # zero-order displacement & first-order velocity overshooting
        # behavior and discontinuous acceleration
        # This is the Wood-Bossak-Zienkiewicz method (Wood, Bossak &
        # Zienkiewicz, 1980)
        # rinf must belong to [0 1]
        if rinf<0:
            rinf=0
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 0')

        if rinf>1:
            rinf=1
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1')

        w1   = -15
        w2   =  45
        w3   = -35
        W1L1 = 1.0
        W2L2 = 1/2
        W3L3 = 1/(1+rinf)**2
        W1L4 = 1.0
        W2L5 = (3-rinf)/2/(1+rinf)
        W1L6 = 2/(1+rinf)
        l1 = 1
        l2 = 1/2
        l3 = 1/(1+rinf)**2
        l4 = 1
        l5=(3-rinf)/2/(1+rinf)

    elif alg == "U1-V0-Opt":
        # first-order displacement & zero-order velocity overshooting
        # behavior and optimal numerical dissipation and dispersion
        # rinf must belong to [0 1]
        if rinf<0:
            rinf=0
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 0')

        if rinf>1:
            rinf=1; # mid-point rule a-form algorithm
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1')

        w1   = -30*(3-8*rinf+6*rinf**2)/(9-22*rinf+19*rinf**2)
        w2   =  15*(25-74*rinf+53*rinf**2)/2/(9-22*rinf+19*rinf**2)
        w3   = -35*(3-10*rinf+7*rinf**2)/(9-22*rinf+19*rinf**2)
        W1L1 = (3-rinf)/2/(1+rinf)
        W2L2 =  1/(1+rinf)**2
        W3L3 =  1/(1+rinf)**3
        W1L4 = (3-rinf)/2/(1+rinf)
        W2L5 =  2/(1+rinf)**3
        W1L6 = (2-rinf)/(1+rinf)
        l1 = 1.0
        l2 = 1/2
        l3 = 1/2/(1+rinf)
        l4 = 1.0
        l5 = 1/(1+rinf)

    elif alg == "U1-V0-CA":
        # first-order displacement & zero-order velocity overshooting
        # behavior and continuous acceleration
        # rinf must belong to [1/2 1]
        if rinf<1/2:
            rinf=1/2
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1/2')

        if rinf>1:
            rinf=1; # Newmark average acceleration a-form algorithm
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1')

        w1   = -60*(2-8*rinf+7*rinf**2)/(11-48*rinf+41*rinf**2)
        w2   =  15*(37-140*rinf+127*rinf**2)/2/(11-48*rinf+41*rinf**2)
        w3   = -35*(5-18*rinf+17*rinf**2)/(11-48*rinf+41*rinf**2)
        W1L1 = (1+3*rinf)/2/(1+rinf)
        W2L2 = 2*rinf/(1+rinf)**2
        W3L3 = 2*rinf/(1+rinf)**3
        W1L4 = (1+3*rinf)/2/(1+rinf)
        W2L5 = 4*rinf/(1+rinf)**3
        W1L6 = 1.0
        l1   = 1.0
        l2   = 1/2
        l3   = 1/2/(1+rinf)
        l4   = 1.0
        l5   = 1/(1+rinf)
    elif alg == "U1-V0-DA":
        # first-order displacement & zero-order velocity overshooting behavior
        # and discontinuous acceleration
        # This is the Newmark average acceleration a-form algorithm
        # rinf must belong to [0 1]
        if rinf<0:
            rinf=0
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 0')

        if rinf>1:
            rinf=1
            warnings.warn('Minimum absolute eigenvalue of amplification matrix is set to 1')

        w1 = -30*(3-4*rinf)/(9-11*rinf)
        w2 =  15*(25-37*rinf)/2/(9-11*rinf)
        w3 = -35*(3-5*rinf)/(9-11*rinf)
        W1L1 = (3+rinf)/2/(1+rinf)
        W2L2 = 1/(1+rinf)
        W3L3 = 1/(1+rinf)**2
        W1L4 = (3+rinf)/2/(1+rinf)
        W2L5 = 2/(1+rinf)**2
        W1L6 = 2/(1+rinf)
        l1   = 1.0
        l2   = 1/2
        l3   = 1/(1+rinf)**2
        l4   = 1.0
        l5   = (3-rinf)/2/(1+rinf)

    elif alg == "Newmark ACA":
        # Newmark Average Constant Acceleration method
        w1 = -15
        w2 =  45
        w3 = -35.0
        W1L1 = 1
        W2L2 = 0.25
        W3L3 = 0.25
        W1L4 = 0.5
        W2L5 = 0.5
        W1L6 = 1.0
        l1 =  1.0
        l2 = 1/2  #
        l3 = 1/4 # beta
        l4 = 1.0
        l5 = 1/2 # gamma

    elif alg == "Newmark LA":
        # Newmark Linear Acceleration method
        w1 = -15
        w2 =  45
        w3 = -35
        W1L1 = 1.0
        W2L2 = 1/6
        W3L3 = 1/6
        W1L4 = 1/2
        W2L5 = 1/2
        W1L6 = 1.0
        l1 = 1.0
        l2 = 1/2
        l3 = 1/6 # beta
        l4 = 1.0
        l5 = 1/2 # gamma

    elif alg == "Newmark BA":
        # Newmark Backward Acceleration method
        w1 = -15
        w2 =  45
        w3 = -35
        W1L1 = 1
        W2L2 = 0.5
        W3L3 = 0.5
        W1L4 = 0.5
        W2L5 = 0.5
        W1L6 = 1
        l1 = 1
        l2 = 0.5
        l3 = 0.5
        l4 = 1
        l5 = 0.5
    
    elif alg == "Fox-Goodwin":
        # Fox-Goodwin formula
        w1=-15
        w2 = 45
        w3=-35
        W1L1 = 1
        W2L2 = 1/12
        W3L3 = 1/12
        W1L4 = 0.5
        W2L5 = 0.5
        W1L6 = 1
        l1 = 1.0
        l2 = 0.5
        l3 = 1/12 # beta
        l4 = 1.0
        l5 = 1/2  # gamma
    else:
        raise ValueError('Unknown algorithm.')

    W1 = (1/2 + w1/3 + w2/4 + w3/5)/(1 + w1/2 + w2/3 + w3/4)
#   W2 = (1/3 + w1/4 + w2/5 + w3/6)/(1 + w1/2 + w2/3 + w3/4)
#   W3 = (1/4 + w1/5 + w2/6 + w3/7)/(1 + w1/2 + w2/3 + w3/4)

    return W1, W1L1, W2L2, W3L3, W1L4, W2L5, W1L6, l1, l2, l3, l4, l5