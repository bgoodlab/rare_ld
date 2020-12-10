import numpy
import gzip
import pylab
import sys
import ld_theory

# Process data

NRs = numpy.array([10,30,100,300,1000])
fstars = numpy.logspace(-2.5,-1,20)

sigmasquareds = {fstar:[] for fstar in fstars}
for NR in NRs:

    filename = 'coalescent_output_%d.txt.gz' % NR
    file = gzip.GzipFile(filename,"r")

    
    n11s = []
    n10s = []
    n01s = []
    n00s = []

    weights_x = []
    weights_y = []

    for line in file:
        items = line.split()
        nAB = long(items[0])
        nAb = long(items[1])
        naB = long(items[2])
        nab = long(items[3])
        Tx = float(items[4])
        Ty = float(items[5])
    
        n = nAB+nAb+naB+nab
    
        weights_x.append(Tx)
        weights_y.append(Ty)
        n11s.append(nAB)
        n10s.append(nAb)
        n01s.append(naB)
        n00s.append(nab)
        
    n11s = numpy.array(n11s)
    n10s = numpy.array(n10s)
    n01s = numpy.array(n01s)
    n00s = numpy.array(n00s)
    
    weights_x = numpy.array(weights_x)
    weights_y = numpy.array(weights_y)

    double_weights = weights_x*weights_y
    ns = n11s+n10s+n01s+n00s
    
    for fstar in fstars:
    # Now do version based on counts
            # First calculate numerator
            rsquared_numerators = n11s*(n11s-1)*n00s*(n00s-1)*numpy.power(1-2.0/ns/fstar,n11s-2)*numpy.power(1-1.0/ns/fstar,n10s-0+n01s-0)
            
            rsquared_numerators += -2*n10s*n01s*n11s*n00s*numpy.power(1-2.0/ns/fstar,n11s-1)*numpy.power(1-1.0/ns/fstar,n10s-1+n01s-1)
            rsquared_numerators += n10s*(n10s-1)*n01s*(n01s-1)*numpy.power(1-2.0/ns/fstar,n11s-0)*numpy.power(1-1.0/ns/fstar,n10s-2+n01s-2)
    
            # Divide by sample size
            rsquared_numerators = rsquared_numerators*(ns>3.5)*1.0/(ns*(ns-1)*(ns-2)*(ns-3)+10*(ns<3.5))
    
            #1
            rsquared_denominators = n10s*(n10s-1)*n01s*(n01s-1)*numpy.power(1-2.0/ns/fstar,n11s-0)*numpy.power(1-1.0/ns/fstar,n10s-2+n01s-2)
            #2
            rsquared_denominators += n10s*n01s*(n01s-1)*n00s*numpy.power(1-2.0/ns/fstar,n11s-0)*numpy.power(1-1.0/n/fstar,n10s-1+n01s-2)       
            #3
            rsquared_denominators += n10s*(n10s-1)*n01s*n11s*numpy.power(1-2.0/ns/fstar,n11s-1)*numpy.power(1-1.0/ns/fstar,n10s-2+n01s-1)
            #4
            rsquared_denominators += n10s*n01s*n11s*n00s*numpy.power(1-2.0/ns/fstar,n11s-1)*numpy.power(1-1.0/ns/fstar,n10s-1+n01s-1)           
            #5
            rsquared_denominators += n10s*(n10s-1)*n01s*n00s*numpy.power(1-2.0/ns/fstar,n11s-0)*numpy.power(1-1.0/ns/fstar,n10s-2+n01s-1)       
            #6
            rsquared_denominators += n10s*n01s*n00s*(n00s-1)*numpy.power(1-2.0/ns/fstar,n11s-0)*numpy.power(1-1.0/ns/fstar,n10s-1+n01s-1)
            #7
            rsquared_denominators += n10s*(n10s-1)*n11s*n00s*numpy.power(1-2.0/n/fstar,n11s-1)*numpy.power(1-1.0/n/fstar,n10s-2+n01s-0)
            #8
            rsquared_denominators += n10s*n11s*n00s*(n00s-1)*numpy.power(1-2.0/ns/fstar,n11s-1)*numpy.power(1-1.0/ns/fstar,n10s-1+n01s-0)
            #9
            rsquared_denominators += n10s*n01s*(n01s-1)*n11s*numpy.power(1-2.0/ns/fstar,n11s-1)*numpy.power(1-1.0/ns/fstar,n10s-1+n01s-2)
            #10
            rsquared_denominators += n01s*(n01s-1)*n11s*n00s*numpy.power(1-2.0/ns/fstar,n11s-1)*numpy.power(1-1.0/ns/fstar,n10s-0+n01s-2)
            #11
            rsquared_denominators += n10s*n01s*n11s*(n11s-1)*numpy.power(1-2.0/ns/fstar,n11s-2)*numpy.power(1-1.0/ns/fstar,n10s-1+n01s-1)
            #12
            rsquared_denominators += n01s*n11s*(n11s-1)*n00s*numpy.power(1-2.0/ns/fstar,n11s-2)*numpy.power(1-1.0/ns/fstar,n10s-0+n01s-1)
            #13
            rsquared_denominators += n10s*n01s*n11s*n00s*numpy.power(1-2.0/ns/fstar,n11s-1)*numpy.power(1-1.0/ns/fstar,n10s-1+n01s-1)
            #14
            rsquared_denominators += n01s*n11s*n00s*(n00s-1)*numpy.power(1-2.0/n/fstar,n11s-1)*numpy.power(1-1.0/n/fstar,n10s-0+n01s-1)
            #15
            rsquared_denominators += n10s*n11s*(n11s-1)*n00s*numpy.power(1-2.0/n/fstar,n11s-2)*numpy.power(1-1.0/n/fstar,n10s-1+n01s-0)
            #16
            rsquared_denominators += n11s*(n11s-1)*n00s*(n00s-1)*numpy.power(1-2.0/ns/fstar,n11s-2)*numpy.power(1-1.0/ns/fstar,n10s-0+n01s-0)
    
            # divide by sample size
            rsquared_denominators = rsquared_denominators*(ns>3.5)*1.0/(ns*(ns-1)*(ns-2)*(ns-3)+10*(ns<3.5))
    
            sigmasquareds[fstar].append( (rsquared_numerators*double_weights).mean()/(rsquared_denominators*double_weights).mean() )
            
for fstar in fstars:
        rhos = 2*NRs*fstar
        
        collapse_xs = 2*NRs*fstar
        collapse_ys = sigmasquareds[fstar]/fstar
        line, = pylab.loglog(collapse_xs,collapse_ys,'o',markersize=3)
        color = pylab.getp(line,'color')

pylab.figure(1)    
theory_xs = numpy.logspace(-3,2,30)
theory_ys = numpy.array([ld_theory.scaled_neutral_ld(x) for x in theory_xs])
pylab.loglog(theory_xs, theory_ys,'k-',linewidth=2)

pylab.ylim([1e-02,6])
        
pylab.savefig('coalescent_sigmasquared_collapse.pdf',bbox_inches='tight')