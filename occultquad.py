from numpy import size,zeros,where,arccos,sqrt,pi,log

# Computes Hasting's polynomial approximation for the complete
# elliptic integral of the first (ek) and second (kk) kind
def ellke(k):
    m1=1.-k**2
    logm1 = log(m1)

    a1=0.44325141463
    a2=0.06260601220
    a3=0.04757383546
    a4=0.01736506451
    b1=0.24998368310
    b2=0.09200180037
    b3=0.04069697526
    b4=0.00526449639
    ee1=1.+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
    ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*(-logm1)
    ek = ee1+ee2
        
    a0=1.38629436112
    a1=0.09666344259
    a2=0.03590092383
    a3=0.03742563713
    a4=0.01451196212
    b0=0.5
    b1=0.12498593597
    b2=0.06880248576
    b3=0.03328355346
    b4=0.00441787012
    ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
    ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*logm1
    kk = ek1-ek2
    
    return [ek,kk]

# Computes the complete elliptical integral of the third kind using
# the algorithm of Bulirsch (1965):
def ellpic_bulirsch(n,k):
    kc=sqrt(1.-k**2); la=n+1.
    if(min(la) < 0.):
        print('Negative l')
    m0=1.; c=1.; la=sqrt(la); d=1./la; e=kc
    while 1:
        f = c; c = d/la+c; g = e/la; d = 2.*(f*g+d)
        la = g + la; g = m0; m0 = kc + m0
        if max(abs(1.-kc/g)) > 1.e-8:
            kc = 2*sqrt(e); e=kc*m0
        else:
            return 0.5*pi*(c*m0+d)/(m0*(m0+la))

#   Python translation of IDL code.
#   This routine computes the lightcurve for occultation of a
#   quadratically limb-darkened source without microlensing.  Please
#   cite Mandel & Agol (2002) and Eastman & Agol (2008) if you make use
#   of this routine in your research.  Please report errors or bugs to
#   jdeast@astronomy.ohio-state.edu
def occultquad(z,u1,u2,p0):

    nz = size(z)
    lambdad = zeros(nz)
    etad = zeros(nz)
    lambdae = zeros(nz)
    omega=1.-u1/3.-u2/6.

    ## tolerance for double precision equalities
    ## special case integrations
    tol = 1e-14

    p = abs(p0)
    
    z = where(abs(p-z) < tol,p,z)
    z = where(abs((p-1)-z) < tol,p-1.,z)
    z = where(abs((1-p)-z) < tol,1.-p,z)
    z = where(z < tol,0.,z)
               
    x1=(p-z)**2.
    x2=(p+z)**2.
    x3=p**2.-z**2.
    
    ## trivial case of no planet
    if p <= 0.:
        muo1 = zeros(nz) + 1. 
        mu0  = zeros(nz) + 1.
        return [muo1,mu0]

    ## Case 1 - the star is unocculted:
    ## only consider points with z lt 1+p
    notusedyet = where( z < (1. + p) )
    notusedyet = notusedyet[0]
    if size(notusedyet) == 0:
        muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*(p > z))+ \
                  u2*etad)/omega
        mu0=1.-lambdae
        return [muo1,mu0]

    # Case 11 - the  source is completely occulted:
    if p >= 1.:
        occulted = where(z[notusedyet] <= p-1.)#,complement=notused2)
        if size(occulted) != 0:
            ndxuse = notusedyet[occulted]
            etad[ndxuse] = 0.5 # corrected typo in paper
            lambdae[ndxuse] = 1.
            # lambdad = 0 already
            notused2 = where(z[notusedyet] > p-1)
            if size(notused2) == 0:
                muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.* \
                                                 (p > z))+u2*etad)/omega
                mu0=1.-lambdae
                return [muo1,mu0]
            notusedyet = notusedyet[notused2]
                
    # Case 2, 7, 8 - ingress/egress (uniform disk only)
    inegressuni = where((z[notusedyet] >= abs(1.-p)) & (z[notusedyet] < 1.+p))
    if size(inegressuni) != 0:
        ndxuse = notusedyet[inegressuni]
        tmp = (1.-p**2.+z[ndxuse]**2.)/2./z[ndxuse]
        tmp = where(tmp > 1.,1.,tmp)
        tmp = where(tmp < -1.,-1.,tmp)
        kap1 = arccos(tmp)
        tmp = (p**2.+z[ndxuse]**2-1.)/2./p/z[ndxuse]
        tmp = where(tmp > 1.,1.,tmp)
        tmp = where(tmp < -1.,-1.,tmp)
        kap0 = arccos(tmp)
        tmp = 4.*z[ndxuse]**2-(1.+z[ndxuse]**2-p**2)**2
        tmp = where(tmp < 0,0,tmp)
        lambdae[ndxuse] = (p**2*kap0+kap1 - 0.5*sqrt(tmp))/pi
        # eta_1
        etad[ndxuse] = 1./2./pi*(kap1+p**2*(p**2+2.*z[ndxuse]**2)*kap0- \
           (1.+5.*p**2+z[ndxuse]**2)/4.*sqrt((1.-x1[ndxuse])*(x2[ndxuse]-1.)))
    
    # Case 5, 6, 7 - the edge of planet lies at origin of star
    ocltor = where(z[notusedyet] == p)#, complement=notused3)
    t = where(z[notusedyet] == p)
    if size(ocltor) != 0:
        ndxuse = notusedyet[ocltor] 
        if p < 0.5:
            # Case 5
            q=2.*p  # corrected typo in paper (2k -> 2p)
            Ek,Kk = ellke(q)
            # lambda_4
            lambdad[ndxuse] = 1./3.+2./9./pi*(4.*(2.*p**2-1.)*Ek+\
                                              (1.-4.*p**2)*Kk)
            # eta_2
            etad[ndxuse] = p**2/2.*(p**2+2.*z[ndxuse]**2)        
            lambdae[ndxuse] = p**2 # uniform disk
        elif p > 0.5:
            # Case 7
            q=0.5/p # corrected typo in paper (1/2k -> 1/2p)
            Ek,Kk = ellke(q)
            # lambda_3
            lambdad[ndxuse] = 1./3.+16.*p/9./pi*(2.*p**2-1.)*Ek-\
                              (32.*p**4-20.*p**2+3.)/9./pi/p*Kk
            # etad = eta_1 already
        else:
            # Case 6
            lambdad[ndxuse] = 1./3.-4./pi/9.
            etad[ndxuse] = 3./32.
        notused3 = where(z[notusedyet] != p)
        if size(notused3) == 0:
            muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*\
                      (lambdad+2./3.*(p > z))+u2*etad)/omega
            mu0=1.-lambdae
            return [muo1,mu0]
        notusedyet = notusedyet[notused3]

    # Case 2, Case 8 - ingress/egress (with limb darkening)
    inegress = where( ((z[notusedyet] > 0.5+abs(p-0.5)) & \
                       (z[notusedyet] < 1.+p))  | \
                      ( (p > 0.5) & (z[notusedyet] > abs(1.-p)) & \
                        (z[notusedyet] < p)) )#, complement=notused4)
    if size(inegress) != 0:

        ndxuse = notusedyet[inegress]
        q=sqrt((1.-x1[ndxuse])/(x2[ndxuse]-x1[ndxuse]))
        Ek,Kk = ellke(q)
        n=1./x1[ndxuse]-1.

        # lambda_1:
        lambdad[ndxuse]=2./9./pi/sqrt(x2[ndxuse]-x1[ndxuse])*\
                         (((1.-x2[ndxuse])*(2.*x2[ndxuse]+x1[ndxuse]-3.)-\
                           3.*x3[ndxuse]*(x2[ndxuse]-2.))*Kk+(x2[ndxuse]-\
                           x1[ndxuse])*(z[ndxuse]**2+7.*p**2-4.)*Ek-\
                          3.*x3[ndxuse]/x1[ndxuse]*ellpic_bulirsch(n,q))

        notused4 = where( ( (z[notusedyet] <= 0.5+abs(p-0.5)) | \
                            (z[notusedyet] >= 1.+p) ) & ( (p <= 0.5) | \
                            (z[notusedyet] <= abs(1.-p)) | \
                            (z[notusedyet] >= p) ))
        if size(notused4) == 0:
            muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*\
                                                     (p > z))+u2*etad)/omega
            mu0=1.-lambdae
            return [muo1,mu0]
        notusedyet = notusedyet[notused4]

    # Case 3, 4, 9, 10 - planet completely inside star
    if p < 1.:
        inside = where(z[notusedyet] <= (1.-p))#, complement=notused5)
        if size(inside) != 0:
            ndxuse = notusedyet[inside]

            ## eta_2
            etad[ndxuse] = p**2/2.*(p**2+2.*z[ndxuse]**2)

            ## uniform disk
            lambdae[ndxuse] = p**2

            ## Case 4 - edge of planet hits edge of star
            edge = where(z[ndxuse] == 1.-p)#, complement=notused6)
            if size(edge[0]) != 0:
                ## lambda_5
                lambdad[ndxuse[edge]] = 2./3./pi*arccos(1.-2.*p)-\
                                      4./9./pi*sqrt(p*(1.-p))*(3.+2.*p-8.*p**2)
                if p > 0.5:
                    lambdad[ndxuse[edge]] -= 2./3.
                notused6 = where(z[ndxuse] != 1.-p)
                if size(notused6) == 0:
                    muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*\
                              (lambdad+2./3.*(p > z))+u2*etad)/omega
                    mu0=1.-lambdae
                    return [muo1,mu0]
                ndxuse = ndxuse[notused6[0]]

            ## Case 10 - origin of planet hits origin of star
            origin = where(z[ndxuse] == 0)#, complement=notused7)
            if size(origin) != 0:
                ## lambda_6
                lambdad[ndxuse[origin]] = -2./3.*(1.-p**2)**1.5
                notused7 = where(z[ndxuse] != 0)
                if size(notused7) == 0:
                    muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*\
                              (lambdad+2./3.*(p > z))+u2*etad)/omega
                    mu0=1.-lambdae
                    return [muo1,mu0]
                ndxuse = ndxuse[notused7[0]]
   
            q=sqrt((x2[ndxuse]-x1[ndxuse])/(1.-x1[ndxuse]))
            n=x2[ndxuse]/x1[ndxuse]-1.
            Ek,Kk = ellke(q)    

            ## Case 3, Case 9 - anywhere in between
            ## lambda_2
            lambdad[ndxuse] = 2./9./pi/sqrt(1.-x1[ndxuse])*\
                              ((1.-5.*z[ndxuse]**2+p**2+x3[ndxuse]**2)*Kk+\
                               (1.-x1[ndxuse])*(z[ndxuse]**2+7.*p**2-4.)*Ek-\
                               3.*x3[ndxuse]/x1[ndxuse]*ellpic_bulirsch(n,q))

        ## if there are still unused elements, there's a bug in the code
        ## (please report it)
        notused5 = where(z[notusedyet] > (1.-p))
        if notused5[0] != 0:
            print("ERROR: the following values of z didn't fit into a case:")
            return [-1,-1]

        muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*(p > z))+\
                  u2*etad)/omega
        mu0=1.-lambdae
        return [muo1,mu0]
