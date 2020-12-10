import numpy
from scipy.integrate import quad
from math import log,exp
from scipy.special import exprel

def ohta_sigmasquared(rho):

    return (10+rho)/(22+13*rho+rho*rho)

def single_locus_sfs(x, Ns):
    
    if Ns>1e-03:
    
        return 2*Ns*exprel(-2*Ns*(1-x))/x
        
    elif Ns<-1e-03:
    
        absNs = -1*Ns
        return (numpy.exp(-2*absNs*x)-numpy.exp(-2*absNs))/x/(1-x)
        #return exprel(2*absNs*(1-x))/x*2*absNs*exp(-2*absNs)
        
    else:
        return 1.0/x
        
def calculate_theta(N,s):
    
    xs = numpy.arange(1,N)*1.0/N
    sfs = single_locus_sfs(xs,N*s)
    
    return sfs.sum()/N

def unscaled_selection_ld(rho,gammaA,gammaB,eps):
    return (gammaA+gammaB+rho)/numpy.square(gammaA+gammaB+eps+rho)
    
def scaled_asexual_ld(gammaA,gammaB,eps,k=2):

    gammaA*=1.0
    gammaB*=1.0
    eps*=1.0
    gamma=gammaA+gammaB+eps

    if k==2:

        if gammaA<(gamma*0.01) and gammaB<(gamma*0.01):
            integrand = lambda z: numpy.exp(-z)*(1-numpy.exp(-z))/numpy.power(gamma/2.0+1-numpy.exp(-z),3)*( 2*(gamma+z))/4.0
    
        elif gammaA<(gamma*0.01):
        
            integrand = lambda z: numpy.exp(-z)*(1-numpy.exp(-z))/numpy.power(gamma/2.0+1-numpy.exp(-z),3)*( (gammaB+1)*(gamma+z) + gamma*(gammaA+1)*(gammaB+1-numpy.exp(-gammaB/gamma*z))/gammaB)/4.0
    
        elif gammaB<(gamma*0.01):
        
            integrand = lambda z: numpy.exp(-z)*(1-numpy.exp(-z))/numpy.power(gamma/2.0+1-numpy.exp(-z),3)*( (gammaA+1)*(gamma+z) + gamma*(gammaB+1)*(gammaA+1-numpy.exp(-gammaA/gamma*z))/gammaA)/4.0
    
        
        else:
        
            integrand = lambda z: numpy.exp(-z)*(1-numpy.exp(-z))/numpy.power(gamma/2.0+1-numpy.exp(-z),3)*( gamma*(gammaB+1)*(gammaA+1-numpy.exp(-gammaA/gamma*z))/gammaA + gamma*(gammaA+1)*(gammaB+1-numpy.exp(-gammaB/gamma*z))/gammaB)/4.0
    
    elif k==1:
        
        if gammaA<(gamma*0.01) and gammaB<(gamma*0.01):
            integrand = lambda z: -numpy.exp(-z)+numpy.exp(-z)/numpy.power(gamma/2.0+1-numpy.exp(-z),2)*( 2*(gamma+z))/4.0
    
        elif gammaA<(gamma*0.01):
        
            integrand = lambda z: -numpy.exp(-z)+numpy.exp(-z)/numpy.power(gamma/2.0+1-numpy.exp(-z),2)*( (gammaB+1)*(gamma+z) + gamma*(gammaA+1)*(gammaB+1-numpy.exp(-gammaB/gamma*z))/gammaB)/4.0
    
        elif gammaB<(gamma*0.01):
        
            integrand = lambda z: -numpy.exp(-z)+numpy.exp(-z)/numpy.power(gamma/2.0+1-numpy.exp(-z),2)*( (gammaA+1)*(gamma+z) + gamma*(gammaB+1)*(gammaA+1-numpy.exp(-gammaA/gamma*z))/gammaA)/4.0
    
        
        else:
        
            integrand = lambda z: -numpy.exp(-z)+numpy.exp(-z)/numpy.power(gamma/2.0+1-numpy.exp(-z),2)*( gamma*(gammaB+1)*(gammaA+1-numpy.exp(-gammaA/gamma*z))/gammaA + gamma*(gammaA+1)*(gammaB+1-numpy.exp(-gammaB/gamma*z))/gammaB)/4.0
    
    result, error = quad(integrand,0,numpy.inf)
    #print result,error
    
    return result

def scaled_asexual_sigmasquared(gamma):

    integrand = lambda z: numpy.exp(-z)*(1-numpy.exp(-z))/numpy.power(gamma/2.0+1-numpy.exp(-z),3)*( 2*(gamma+z))/4.0
    
    result, error = quad(integrand,0,numpy.inf)
    #print result,error
    return result
    
def scaled_asexual_sigmafourth(gamma):
    integrand = lambda z: numpy.exp(-z)*numpy.power(1-numpy.exp(-z),3)/numpy.power(gamma/2.0+1-numpy.exp(-z),5)*( 2*(gamma+z))/4.0
    
    result, error = quad(integrand,0,numpy.inf)
    #print result,error
    return result

def asexual_eta(gamma):
    
    return scaled_asexual_sigmafourth(gamma)/scaled_asexual_sigmasquared(gamma)
        
from scipy.special import expi as Ei

def scaled_neutral_le(rho):
    
    if rho>100:
        return 1.0
    
    
    A = lambda z: rho*numpy.exp(rho)*(Ei(-rho-z)-Ei(-rho))
    N = lambda z: numpy.exp(-z)
    D1 = lambda z: (rho+z)*(2+rho+z)
    D = lambda z: (rho+z)*(2+rho+z)*(-1+A(z))+rho*(1+rho+z)*numpy.exp(-z)
    
    integrand = lambda z: rho*N(z)/numpy.square(D(z))
    
    result, error = quad(integrand,0,numpy.inf, limit=100)
    
    #print result,error
    return result

def scaled_neutral_small_ld(rho):
    
    if rho>100:
        return 1/numpy.power(rho,1)
    
    A = lambda z: rho*numpy.exp(rho)*(Ei(-rho-z)-Ei(-rho))
    N = lambda z: 2*numpy.exp(-z)
    D1 = lambda z: (rho+z)*(2+rho+z)
    D = lambda z: (rho+z)*(2+rho+z)*(-rho+1.0+rho*A(z))+rho*rho*(1+rho+z)*numpy.exp(-z)
    
    integrand = lambda z: D1(z)*N(z)/numpy.square(D(z))*(-4*(2+rho)+8*D1(z)/D(z))/rho
         
    result, error = quad(integrand,0,numpy.inf, limit=100)
    
    #print result,error
    return result

def scaled_neutral_ld(rho,k=2):
    
    if rho>100:
        if k==2:
            return 1/numpy.power(rho,1)
        else:
            return 2/numpy.power(rho,3)
    
    
    A = lambda z: rho*numpy.exp(rho)*(Ei(-rho-z)-Ei(-rho))
    N = lambda z: numpy.exp(-z)
    D1 = lambda z: (rho+z)*(2+rho+z)
    D = lambda z: (rho+z)*(2+rho+z)*(-1+A(z))+rho*(1+rho+z)*numpy.exp(-z)
    if k==2:
        integrand = lambda z: D1(z)*N(z)/numpy.square(D(z))*((2+rho)+D1(z)/D(z))
        
    elif k==4:
        
        integrand = lambda z: D1(z)*N(z)/numpy.square(D(z))*(numpy.power(2+rho,3)+3*D1(z)*numpy.square(2+rho)/D(z)+3*numpy.square(D1(z))*(2+rho)/numpy.square(D(z))+numpy.power(D1(z),3)/numpy.power(D(z),3))
    
    elif k==1:
        integrand = lambda z: -N(z)+D1(z)*N(z)/numpy.square(D(z))
    else:
        integrand = lambda z: 0
            
    result, error = quad(integrand,0,numpy.inf, limit=100)
    
    #print result,error
    return result

def approx_scaled_neutral_ld(rho,k=2):
    
    if rho>10:
        return 1
        
    
    A = lambda z: rho*numpy.exp(rho)*(Ei(-rho-z)-Ei(-rho))
    N = lambda z: numpy.exp(-z)
    D1 = lambda z: (rho+z)*(2+z)
    D = lambda z: -1*(rho+z)*(2+z)+rho*(1+z)*numpy.exp(-z)
    
    if k==2:
        integrand = lambda z: D1(z)*N(z)/numpy.square(D(z))*((2)+D1(z)/D(z))
        
    elif k==4:
        
        integrand = lambda z: D1(z)*N(z)/numpy.square(D(z))*(numpy.power(2,3)+3*D1(z)*numpy.square(2)/D(z)+3*numpy.square(D1(z))*(2)/numpy.square(D(z))+numpy.power(D1(z),3)/numpy.power(D(z),3))
    
    else:
        integrand = lambda z: 0
            
    result, error = quad(integrand,0,numpy.inf, limit=100)
    
    #print result,error
    return result
    
    
def neutral_eta(rho):
    
    return scaled_neutral_ld(rho,k=4)/scaled_neutral_ld(rho,k=2)

def qle_neutral_eta(rho,fstar):

    weight = numpy.exp(-1.0/rho/fstar)
    return neutral_eta(rho)*(1-weight)+(2+rho*fstar)/rho/rho*weight
    
    
def approx_neutral_eta(rho):
    
    return approx_scaled_neutral_ld(rho,k=4)/approx_scaled_neutral_ld(rho,k=2)


def approx_eta(rho):
    
    return 1.0/(1+(rho**2/2.0))

def very_approx_eta(rho):
    
    return (1.0/3)*numpy.log(1+1.0/numpy.power(rho/(6.0*4)**0.5,3))/numpy.log(1+4*(6.0*4)**0.5/rho)
    
    # a^3*c/b = 2
    # 3c = 1-->c=1/3 --> a^3/b = 6 --> a^3 = 6 b --> a = (6 b)^(1/3)
    #
    # a = b/2 --> (a)^2 = 12 --> a = sqrt(6*const), b=const*a
    #
    # but also need a<b --> a=b --> a^2 = 6 --> a = sqrt(6), b=sqrt(6)
    
def approx_ld(rho):
    return 1.0/2.0/((1+rho/2.0)**2)*((1+rho/4.0)*log(1+2.0/rho)+0.5)*(1+rho/2)
    
def very_approx_ld(rho):
    return 0.5*log(1+2.0/rho)

def bad_ld(rho):
    return 5.0/13*log(1+13.0/5/rho)
    
if __name__=='__main__':

    import pylab
    import sys
    
    rhos = numpy.logspace(-10,2.1,50)
    lds = numpy.array([scaled_neutral_ld(rho) for rho in rhos])
    approx_lds = numpy.array([0.5*log(2.0/rho) for rho in rhos])
    pylab.loglog(rhos,lds,'-')
    pylab.loglog(rhos,approx_lds,':')
    pylab.show()
    sys.exit(0)
    
    N = 1e05
    print calculate_theta(N,0), log(N)
    print calculate_theta(N,-1e-02)
    
    pylab.figure()
    
    gammas = numpy.logspace(-3,2,50)
    
    theory_ys = numpy.array([scaled_asexual_ld(0,0,gamma) for gamma in gammas])
    
    
    pylab.loglog(gammas,theory_ys,'k-')
    
    theory_ys = numpy.array([scaled_asexual_ld(gamma/2,gamma/2,0) for gamma in gammas])
    
    
    pylab.loglog(gammas,theory_ys,'r-')
    
    theory_ys = numpy.array([scaled_asexual_ld(gamma,0,0) for gamma in gammas])
    pylab.loglog(gammas,theory_ys,'b-')
    
    theory_ys = numpy.array([ohta_sigmasquared(gamma) for gamma in gammas])
    pylab.semilogx(gammas,theory_ys,'b-')
    
    pylab.semilogx(gammas,numpy.zeros_like(gammas),'k:')
    
    
    pylab.ylim([1e-02,10])
    
    
    pylab.figure(2)
    gammas = numpy.logspace(-3,2,50)
    
    theory_ys = numpy.array([scaled_asexual_ld(0,0,gamma,k=1) for gamma in gammas])
    
    
    pylab.semilogx(gammas,theory_ys,'k-')
    
    theory_ys = numpy.array([scaled_asexual_ld(gamma/2,gamma/2,0,k=1) for gamma in gammas])
    
    
    pylab.semilogx(gammas,theory_ys,'r-')
    
    theory_ys = numpy.array([scaled_asexual_ld(gamma,0,0,k=1) for gamma in gammas])
    pylab.semilogx(gammas,theory_ys,'b-')
    
    pylab.semilogx(gammas,numpy.zeros_like(gammas),'k:')
    
    # Try to do a bunch of Rs with Ns = 100 so Ns*f = 1e-03
    
    
    
    pylab.show()