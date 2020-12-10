import sys
import numpy
import pylab
import gzip
import parameters
from math import log10
import ld_theory
from math import log

type = "r"

fstars = numpy.array([1e-03,3e-03,1e-02,3e-02,1e-01,1]) #numpy.logspace(-3,0,50)
params = parameters.params[type]

asexual_sigmasquareds = []
sigmasquareds = {fstar:[] for fstar in fstars}
sigmafourths = {fstar:[] for fstar in fstars}

denominatorsquareds = {fstar:[] for fstar in fstars}
denominatorfourths = {fstar:[] for fstar in fstars}

bare_denominatorsquareds = {fstar:[] for fstar in fstars}
bare_denominatorfourths = {fstar:[] for fstar in fstars}


Nrs = []

#for param_idx in xrange(0,len(params)):
for param_idx in [0,5,10]:  
    N = params[param_idx][2]
    r = params[param_idx][6]
    Nr = N*r
    Nmu = 1.0/log(N)
    
    if Nr>0 and Nr<1:
        continue
        
    Nrs.append(N*r)
    
    filename = 'output/output_%s_%d.txt.gz' % (type,param_idx)
    file = gzip.GzipFile(filename,"r")
    f11s = []
    f10s = []
    f01s = []

    print "Loading %s" % filename
    for line in file:
    
        if line.startswith('//'):
            continue

        items = line.split()
        f11s.append(float(items[0]))
        f10s.append(float(items[1]))
        f01s.append(float(items[2]))
    file.close()
   
    f11s = numpy.array(f11s)
    f10s = numpy.array(f10s)
    f01s = numpy.array(f01s)
    f00s = 1-f11s-f10s-f01s

    fAs = f11s+f10s
    fBs = f11s+f01s

    #Ds = f11s
    #sampling_variances = f01s*f10s
    
    Ds = f11s*f00s-f10s*f01s
    sampling_variances = fAs*(1-fAs)*fBs*(1-fBs)

    Dsquareds = numpy.square(Ds)
    
    for fstar in fstars:
        
        Hs = numpy.exp(-fAs/fstar-fBs/fstar)
    
        sigmasquared_numerator = (Dsquareds*Hs).mean() 
        sigmasquared_denominator =  (sampling_variances*Hs).mean() 
        sigmasquared = sigmasquared_numerator/sigmasquared_denominator
        sigmasquareds[fstar].append(sigmasquared)
        
        sigmafourth_numerator = (numpy.square(Dsquareds)*Hs).mean() 
        sigmafourth_denominator =  (numpy.square(sampling_variances)*Hs).mean() 
        sigmafourth = sigmafourth_numerator/sigmafourth_denominator/3.0
        sigmafourths[fstar].append(sigmafourth)
        
        denominatorsquareds[fstar].append((sampling_variances*Hs).mean())
        denominatorfourths[fstar].append((numpy.square(sampling_variances)*Hs).mean())
        
        bare_denominatorsquareds[fstar].append((f10s*f01s*Hs).mean())
        bare_denominatorfourths[fstar].append((numpy.square(f10s*f01s)*Hs).mean())
        
        
        if r==0:
            asexual_sigmasquareds.append(sigmasquared)

for fstar in fstars:
    sigmasquareds[fstar] = numpy.array(sigmasquareds[fstar])
    sigmafourths[fstar] = numpy.array(sigmafourths[fstar])
    
    denominatorsquareds[fstar] = numpy.array(denominatorsquareds[fstar])
    denominatorfourths[fstar] = numpy.array(denominatorfourths[fstar])
    
    bare_denominatorsquareds[fstar] = numpy.array(bare_denominatorsquareds[fstar])
    
    bare_denominatorfourths[fstar] = numpy.array(bare_denominatorfourths[fstar])
    
    
print Nrs
#Nrs[0]=1e-02
Nrs = numpy.array(Nrs)
Nrs = numpy.clip(Nrs,1e-02,1e09)

asexual_sigmasquareds = numpy.array(asexual_sigmasquareds)

for fstar in fstars:

    theory_sigmasquareds = fstar/2.0/numpy.square(1+Nrs*fstar)*((1+Nrs*fstar/2)*numpy.log(1+1.0/numpy.clip(Nrs,1,1e09)/fstar)+0.5)*(1+Nrs*fstar)
    pylab.figure(1)
    line, = pylab.loglog(Nrs,sigmasquareds[fstar],'.',label=('%g' % log10(fstar)))
    color = pylab.getp(line,'color')
    pylab.loglog(Nrs,theory_sigmasquareds,'-',color=color)
    pylab.figure(2)
    collapse_xs = 2*Nrs*fstar
    collapse_ys = sigmasquareds[fstar]/fstar
    good_idxs = (Nrs>=1)
    pylab.loglog(collapse_xs[good_idxs],collapse_ys[good_idxs],'.',color=color)
    pylab.figure(3)
    collapse_ys = sigmasquareds[fstar]*Nrs
    pylab.loglog(collapse_xs[good_idxs],collapse_ys[good_idxs],'.',color=color)
    pylab.figure(4)
    collapse_ys = sigmafourths[fstar]/sigmasquareds[fstar]
    pylab.loglog(collapse_xs[good_idxs],collapse_ys[good_idxs],'.',color=color)
    
    pylab.figure(5)
    
    pylab.loglog(collapse_xs, denominatorsquareds[fstar]/fstar**2/Nmu**2,'o',color=color,alpha=0.5)
    pylab.loglog(collapse_xs, bare_denominatorsquareds[fstar]/fstar**2/Nmu**2,'s',color=color,alpha=0.5)
    
    pylab.loglog(collapse_xs, denominatorfourths[fstar]/fstar**4/Nmu**2,'^',color=color,alpha=0.5)
    pylab.loglog(collapse_xs, bare_denominatorfourths[fstar]/fstar**4/Nmu**2,'v',color=color,alpha=0.5)
    
    
    
pylab.figure(1)    
pylab.legend(loc='lower left',frameon=False)

#frequency_bins = numpy.linspace(0,1,1000)        
#pAs = numpy.histogram(fAs,bins=frequency_bins)[0]
#pBs = numpy.histogram(fBs,bins=frequency_bins)[0]
   
#print len(pAs), len(frequency_bins)  
#pylab.loglog(frequency_bins[1:], pAs,'-')
#pylab.loglog(frequency_bins[1:], pBs,'-')


pylab.savefig('twolocus_recombination_output.pdf',bbox_inches='tight')

pylab.figure(2)

theory_xs = numpy.logspace(-3,4,50)
theory_ys = numpy.array([ld_theory.scaled_neutral_ld(x) for x in theory_xs])
approx_theory_ys = numpy.array([ld_theory.bad_ld(x) for x in theory_xs])
asexual_ys = numpy.array([ld_theory.scaled_asexual_sigmasquared(x) for x in theory_xs])

pylab.loglog(theory_xs, theory_ys,'k-',linewidth=2)
pylab.loglog(theory_xs, approx_theory_ys,'k:',linewidth=2)

theory_xs = numpy.logspace(-3,2,50)
asexual_ys = numpy.array([ld_theory.scaled_asexual_sigmasquared(x) for x in theory_xs])
pylab.loglog(theory_xs,asexual_ys,'r-')

pylab.savefig('twolocus_recombination_collapse_1.pdf',bbox_inches='tight')

pylab.figure(3)
pylab.savefig('twolocus_recombination_collapse_2.pdf',bbox_inches='tight')

pylab.figure(4)

rhos = numpy.logspace(-3,2,40)
theory_etas = numpy.array([ld_theory.neutral_eta(rho) for rho in rhos])
pylab.loglog(rhos,theory_etas,'k-',linewidth=2)   

pylab.savefig('twolocus_recombination_collapse_3.pdf',bbox_inches='tight')

pylab.figure(5)

rhos = numpy.logspace(-3,2,40)
pylab.loglog(rhos,numpy.ones_like(rhos),'k-',linewidth=2)
pylab.legend(loc='lower left',frameon=False)
pylab.savefig('twolocus_recombination_collapse_4.pdf',bbox_inches='tight')

