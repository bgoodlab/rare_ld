import sys
import numpy
import pylab
import gzip
import parameters
from math import log10,log
import ld_theory
from numpy.random import multinomial

import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from numpy.random import randint, shuffle, poisson, binomial, choice, hypergeometric
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

theory_xs = numpy.logspace(-1,3,50)

# Set up figure
f, axes = pylab.subplots(2,1,sharex=True,squeeze=False)
f.set_size_inches(3.42,5)

k1_axis = axes[0][0]
k1_axis.set_ylabel("$\\sigma_d^1$")
k1_axis.semilogx(theory_xs,numpy.zeros_like(theory_xs),'k:')

k2_axis = axes[1][0]
k2_axis.set_ylabel("$\\sigma_d^2$")
k2_axis.set_xlabel("$2 N R$")

# Plot DATA

n=1e05 # actual population size
params = parameters.params

selection_gammaA = 2*params['r_selA'][0][2]*params['r_selA'][0][3]
selection_gammaB = 2*params['r_selA'][0][2]*params['r_selA'][0][4]
selection_eps = 0
selection_fstar = 1.0/selection_gammaA

print "Effective f0 =", selection_fstar

for type,symbol,counts_symbol,color in zip(['r','r_selA'],['o','s'],['.','.'],['k','r']):
    
    
    if type=='r':
        fstars = numpy.array([selection_fstar,100])
    elif type=='r_selA':
        fstars = numpy.array([100])
    else:
        print "type=", type
    
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
            #continue
        
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
            
            if False:
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
            
            else:
                pass
                #sigmasquareds_counts[fstar].append(0)
            
            bare_numeratorsquareds[fstar].append(( numpy.square(bare_Ds)*Hs).mean()/theta**2)

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
        collapse_xs = gammas
        collapse_ys = sigmasquareds[fstar]
        
        if fstar<1:
            colorVal='0.7'
        else:
            colorVal=color  
        
        line, = k2_axis.loglog(collapse_xs,collapse_ys,symbol,markersize=3,color=colorVal)
        
        
        
        #collapse_ys = sigmasquareds_counts[fstar]/fstar
        #pylab.loglog(collapse_xs,collapse_ys,counts_symbol,markersize=3)
        
        collapse_ys = sigmas[fstar]
        k1_axis.semilogx(collapse_xs, collapse_ys,symbol, color=colorVal,markersize=3)
        
        
        #pylab.figure(2)
        #collapse_ys = sigmafourths[fstar]/sigmasquareds[fstar]
        #pylab.loglog(collapse_xs,collapse_ys,symbol,color=color,markersize=3)
        
        #pylab.figure(3)
        #collapse_ys = denominatorsquareds[fstar]/fstar**2
        #pylab.semilogx(collapse_xs,collapse_ys,symbol,color=color,markersize=3)
        
        #pylab.figure(4)
        #collapse_ys = bare_numeratorsquareds[fstar]/fstar**3
        #pylab.loglog(collapse_xs,collapse_ys,symbol,color=color,markersize=3)
        
        
theory_ys = numpy.array([ld_theory.ohta_sigmasquared(x) for x in theory_xs])
k2_axis.semilogx(theory_xs,theory_ys,'k-')

theory_ys = numpy.array([ld_theory.scaled_neutral_ld(max([x,1])*selection_fstar)*selection_fstar for x in theory_xs])
k2_axis.loglog(theory_xs, theory_ys,'-',color='0.7')

theory_ys = numpy.array([ld_theory.unscaled_selection_ld(x,selection_gammaA,selection_gammaB,selection_eps) for x in theory_xs])
k2_axis.loglog(theory_xs, theory_ys,'-',color='r')

theory_ys = numpy.array([ld_theory.scaled_neutral_ld(max([x,1])*selection_fstar,k=1) for x in theory_xs])
k1_axis.semilogx(theory_xs, theory_ys,'-',color='0.7')


k2_axis.set_ylim([1e-03,1])
k2_axis.set_xlim([1e-01,1e03])

k1_axis.set_xlim([1e-01,1e03])
k1_axis.set_ylim([-1,2])

k1_axis.plot([1e-05],[1],'o',color='k',label='Neutral, $f_0=\infty$',markersize=3)
k1_axis.plot([1e-05],[1],'o',color='r',label='$2Ns=25$, $f_0=\infty$',markersize=3)
k1_axis.plot([1e-05],[1],'o',color='0.7',label='Neutral, $f_0=1/2Ns$',markersize=3)

k1_axis.legend(frameon=False,loc='upper right',numpoints=1,scatterpoints=1)

pylab.savefig('neutral_selection_comparison.pdf',bbox_inches='tight')
