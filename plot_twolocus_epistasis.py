import sys
import numpy
import pylab
import gzip
import parameters
from math import log10,log
import ld_theory
from numpy.random import multinomial

n=1e05 # actual population size
fstars = numpy.logspace(-2.5,-1,20)
params = parameters.params
for type,symbol,counts_symbol in zip(['eps','r'],['v','^'],['<','>']):
    
    sigmas = {fstar:[] for fstar in fstars}
    sigmasquareds = {fstar:[] for fstar in fstars}
    sigmafourths = {fstar:[] for fstar in fstars}
    denominatorsquareds = {fstar:[] for fstar in fstars}
    bare_numeratorsquareds = {fstar:[] for fstar in fstars}
    
    sigmasquareds_counts = {fstar:[] for fstar in fstars}
    
    gammas = []
    
    for param_idx in xrange(0,len(params[type])):
    
        N = params[type][param_idx][2]
        eps = params[type][param_idx][5]
        r = params[type][param_idx][6]
        gamma = 2*(N*eps+N*r)
        theta = 1.0/log(N)
    
        if gamma<5:
            #print gamma
            pass
            continue
        
        gammas.append(gamma)
    
        filename = 'output/output_%s_%d.txt.gz' % (type,param_idx)
        file = gzip.GzipFile(filename,"r")
        f11s = []
        f10s = []
        f01s = []
        
        # for counts (regenerated on the fly)
        n11s = []
        n10s = []
        n01s = []
        print "Loading %s" % filename
        for line in file:
    
            if line.startswith('//'):
                continue
            items = line.split()
            
            f11 = float(items[0])
            f10 = float(items[1])
            f01 = float(items[2])
            
            f11s.append(f11)
            f10s.append(f10)
            f01s.append(f01)
            
        
            n11 = long(f11*n)
            n10 = long(f10*n)
            n01 = long(f01*n)
            
            n11s.append(n11)
            n10s.append(n10)
            n01s.append(n01)
            
        file.close()
   
        f11s = numpy.array(f11s)
        f10s = numpy.array(f10s)
        f01s = numpy.array(f01s)
        f00s = 1-f11s-f10s-f01s
        fAs = f11s+f10s
        fBs = f11s+f01s

        n11s = numpy.array(n11s)
        n10s = numpy.array(n10s)
        n01s = numpy.array(n01s)
        ns = numpy.ones_like(n11s)*n
        n00s = ns-n10s-n11s-n01s
        
        for fstar in fstars:
        
            bare_Ds = f11s
            bare_variances = f01s*f10s
    
            Ds = f11s*f00s-f10s*f01s
            sampling_variances = fAs*(1-fAs)*fBs*(1-fBs)

            Dsquareds = numpy.square(Ds)
        
            Hs = numpy.exp(-fAs/fstar-fBs/fstar)
    
            sigma_numerator = (Ds*Hs).mean()
            sigma_denominator = (sampling_variances*Hs).mean()
            sigma = sigma_numerator/sigma_denominator
            sigmas[fstar].append(sigma)
    
            sigmasquared_numerator = (Dsquareds*Hs).mean() 
            sigmasquared_denominator =  (sampling_variances*Hs).mean() 
            sigmasquared = sigmasquared_numerator/sigmasquared_denominator
            sigmasquareds[fstar].append(sigmasquared)
        
            sigmafourth_numerator = (numpy.square(Dsquareds)*Hs).mean() 
            sigmafourth_denominator =  (numpy.square(sampling_variances)*Hs).mean() 
            sigmafourth = sigmafourth_numerator/sigmafourth_denominator/3.0
            sigmafourths[fstar].append(sigmafourth)
            
            denominatorsquareds[fstar].append(sigmasquared_denominator/theta**2)
            
            # Now do version based on counts
            # First calculate numerator
            rsquared_numerators = n11s*(n11s-1)*n00s*(n00s-1)*numpy.power(1-2.0/n/fstar,n11s-2)*numpy.power(1-1.0/n/fstar,n10s-0+n01s-0)
            
            rsquared_numerators += -2*n10s*n01s*n11s*n00s*numpy.power(1-2.0/n/fstar,n11s-1)*numpy.power(1-1.0/n/fstar,n10s-1+n01s-1)
            rsquared_numerators += n10s*(n10s-1)*n01s*(n01s-1)*numpy.power(1-2.0/n/fstar,n11s-0)*numpy.power(1-1.0/n/fstar,n10s-2+n01s-2)
    
            # Divide by sample size
            rsquared_numerators = rsquared_numerators*(ns>3.5)*1.0/(ns*(ns-1)*(ns-2)*(ns-3)+10*(ns<3.5))
    
            #1
            rsquared_denominators = n10s*(n10s-1)*n01s*(n01s-1)*numpy.power(1-2.0/n/fstar,n11s-0)*numpy.power(1-1.0/n/fstar,n10s-2+n01s-2)
            #2
            rsquared_denominators += n10s*n01s*(n01s-1)*n00s*numpy.power(1-2.0/n/fstar,n11s-0)*numpy.power(1-1.0/n/fstar,n10s-1+n01s-2)       
            #3
            rsquared_denominators += n10s*(n10s-1)*n01s*n11s*numpy.power(1-2.0/n/fstar,n11s-1)*numpy.power(1-1.0/n/fstar,n10s-2+n01s-1)
            #4
            rsquared_denominators += n10s*n01s*n11s*n00s*numpy.power(1-2.0/n/fstar,n11s-1)*numpy.power(1-1.0/n/fstar,n10s-1+n01s-1)           
            #5
            rsquared_denominators += n10s*(n10s-1)*n01s*n00s*numpy.power(1-2.0/n/fstar,n11s-0)*numpy.power(1-1.0/n/fstar,n10s-2+n01s-1)       
            #6
            rsquared_denominators += n10s*n01s*n00s*(n00s-1)*numpy.power(1-2.0/n/fstar,n11s-0)*numpy.power(1-1.0/n/fstar,n10s-1+n01s-1)
            #7
            rsquared_denominators += n10s*(n10s-1)*n11s*n00s*numpy.power(1-2.0/n/fstar,n11s-1)*numpy.power(1-1.0/n/fstar,n10s-2+n01s-0)
            #8
            rsquared_denominators += n10s*n11s*n00s*(n00s-1)*numpy.power(1-2.0/n/fstar,n11s-1)*numpy.power(1-1.0/n/fstar,n10s-1+n01s-0)
            #9
            rsquared_denominators += n10s*n01s*(n01s-1)*n11s*numpy.power(1-2.0/n/fstar,n11s-1)*numpy.power(1-1.0/n/fstar,n10s-1+n01s-2)
            #10
            rsquared_denominators += n01s*(n01s-1)*n11s*n00s*numpy.power(1-2.0/n/fstar,n11s-1)*numpy.power(1-1.0/n/fstar,n10s-0+n01s-2)
            #11
            rsquared_denominators += n10s*n01s*n11s*(n11s-1)*numpy.power(1-2.0/n/fstar,n11s-2)*numpy.power(1-1.0/n/fstar,n10s-1+n01s-1)
            #12
            rsquared_denominators += n01s*n11s*(n11s-1)*n00s*numpy.power(1-2.0/n/fstar,n11s-2)*numpy.power(1-1.0/n/fstar,n10s-0+n01s-1)
            #13
            rsquared_denominators += n10s*n01s*n11s*n00s*numpy.power(1-2.0/n/fstar,n11s-1)*numpy.power(1-1.0/n/fstar,n10s-1+n01s-1)
            #14
            rsquared_denominators += n01s*n11s*n00s*(n00s-1)*numpy.power(1-2.0/n/fstar,n11s-1)*numpy.power(1-1.0/n/fstar,n10s-0+n01s-1)
            #15
            rsquared_denominators += n10s*n11s*(n11s-1)*n00s*numpy.power(1-2.0/n/fstar,n11s-2)*numpy.power(1-1.0/n/fstar,n10s-1+n01s-0)
            #16
            rsquared_denominators += n11s*(n11s-1)*n00s*(n00s-1)*numpy.power(1-2.0/n/fstar,n11s-2)*numpy.power(1-1.0/n/fstar,n10s-0+n01s-0)
    
            # divide by sample size
            rsquared_denominators = rsquared_denominators*(ns>3.5)*1.0/(ns*(ns-1)*(ns-2)*(ns-3)+10*(ns<3.5))
    
            sigmasquareds_counts[fstar].append( rsquared_numerators.mean()/rsquared_denominators.mean() )
            bare_numeratorsquareds[fstar].append((numpy.square(bare_Ds)*Hs).mean()/theta**2)

    gammas = numpy.array(gammas)
    for fstar in fstars:
        sigmas[fstar] = numpy.array(sigmas[fstar])
        sigmasquareds[fstar] = numpy.array(sigmasquareds[fstar])
        sigmafourths[fstar] = numpy.array(sigmafourths[fstar])
        denominatorsquareds[fstar] = numpy.array(denominatorsquareds[fstar])
        bare_numeratorsquareds[fstar] = numpy.array(bare_numeratorsquareds[fstar])
        
        sigmasquareds_counts[fstar] = numpy.array(sigmasquareds_counts[fstar])
        
        
    for fstar in fstars:

        pylab.figure(1)
        collapse_xs = gammas*fstar
        collapse_ys = sigmasquareds[fstar]/fstar
        line, = pylab.loglog(collapse_xs,collapse_ys,symbol,markersize=3)
        color = pylab.getp(line,'color')
        
        collapse_ys = sigmasquareds_counts[fstar]/fstar
        pylab.loglog(collapse_xs,collapse_ys,counts_symbol,markersize=3)
        
        pylab.figure(2)
        collapse_ys = sigmafourths[fstar]/sigmasquareds[fstar]
        pylab.loglog(collapse_xs,collapse_ys,symbol,color=color,markersize=3)
        
        pylab.figure(3)
        collapse_ys = denominatorsquareds[fstar]/fstar**2
        pylab.semilogx(collapse_xs,collapse_ys,symbol,color=color,markersize=3)
        
        pylab.figure(4)
        collapse_ys = bare_numeratorsquareds[fstar]/fstar**3
        pylab.loglog(collapse_xs,collapse_ys,symbol,color=color,markersize=3)
        
        pylab.figure(5)
        collapse_ys = sigmas[fstar]
        pylab.semilogx(collapse_xs, collapse_ys,symbol, color=color,markersize=3)
        
    
pylab.figure(1)    
theory_xs = numpy.logspace(-3,2,50)
theory_ys = numpy.array([ld_theory.scaled_asexual_sigmasquared(x) for x in theory_xs])
pylab.loglog(theory_xs, theory_ys,'k-',linewidth=2)
theory_xs = numpy.logspace(-3,4,30)
theory_ys = numpy.array([ld_theory.scaled_neutral_ld(x) for x in theory_xs])
pylab.loglog(theory_xs, theory_ys,'k-',linewidth=2)

pylab.savefig('twolocus_epistasis_collapse_sigma.pdf',bbox_inches='tight')

pylab.figure(2)
theory_ys = numpy.array([ld_theory.asexual_eta(x) for x in theory_xs])
pylab.loglog(theory_xs, theory_ys,'k-',linewidth=2)
theory_ys = numpy.array([ld_theory.neutral_eta(x) for x in theory_xs])
pylab.loglog(theory_xs, theory_ys,'k-',linewidth=2)
pylab.savefig('twolocus_epistasis_collapse_eta.pdf',bbox_inches='tight')

pylab.figure(3)    
theory_xs = numpy.logspace(-3,2,50)
theory_ys = numpy.ones_like(theory_xs)
pylab.semilogx(theory_xs, theory_ys,'k-',linewidth=2)
pylab.ylim([0,1.1])
pylab.savefig('twolocus_epistasis_collapse_denominator.pdf',bbox_inches='tight')

pylab.figure(4)    
theory_xs = numpy.logspace(-3,2,50)
theory_ys = numpy.array([ld_theory.scaled_asexual_sigmasquared(x) for x in theory_xs])
pylab.loglog(theory_xs, theory_ys,'k-',linewidth=2)
theory_xs = numpy.logspace(-3,4,50)
theory_ys = numpy.array([ld_theory.scaled_neutral_ld(x) for x in theory_xs])
pylab.loglog(theory_xs, theory_ys,'k-',linewidth=2)

pylab.savefig('twolocus_epistasis_collapse_numerator.pdf',bbox_inches='tight')

pylab.figure(5)
theory_xs = numpy.logspace(-3,2,50)
theory_ys = numpy.array([ld_theory.scaled_neutral_ld(x,k=1) for x in theory_xs])
pylab.semilogx(theory_xs, theory_ys,'k-',linewidth=2)
pylab.savefig('twolocus_epistasis_signed_collapse.pdf',bbox_inches='tight')