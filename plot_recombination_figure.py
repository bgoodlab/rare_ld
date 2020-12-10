import sys
import numpy
import pylab
import gzip
import parameters
from math import log10,log
import ld_theory
from numpy.random import multinomial
from scipy.special import gammaln

import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from numpy.random import randint, shuffle, poisson, binomial, choice, hypergeometric
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

debug = True

mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

theory_xs = numpy.logspace(-3,2.5,50)

# Set up figure
f, axes = pylab.subplots(2,1,sharex=True,squeeze=False)
f.set_size_inches(3.42,5)

k1_axis = axes[0][0]
k1_axis.set_ylabel("$\\sigma_d^1$")
k1_axis.semilogx(theory_xs,numpy.zeros_like(theory_xs),'k:')
k1_axis.set_xlim([theory_xs[0],theory_xs[-1]])

k2_axis = axes[1][0]
k2_axis.set_ylabel("$\\sigma_d^2 / f_0$")
k2_axis.set_xlabel("$2 N R f_0$")
k2_axis.set_xlim([theory_xs[0],theory_xs[-1]])

vmin=-3
vmax=-1
cmap='jet_r'

jet = cm = pylab.get_cmap(cmap) 
cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

pylab.figure(2) # new fig for LE!
le_axis = pylab.gca()
le_axis.set_xlabel("$2 N R f_0$")
le_axis.set_ylabel("Linkage equilibrium, $E(f_0)$")
le_axis.set_xlim([theory_xs[0],theory_xs[-1]])



pylab.figure(3)
denom_axis = pylab.gca()
denom_axis.set_xlabel("$2 N R f_0$")
denom_axis.set_ylabel("Denominator")
denom_axis.set_xlim([theory_xs[0],theory_xs[-1]])

pylab.figure(4) # new fig for LE!
small_ld_axis = pylab.gca()
small_ld_axis.set_xlabel("$2 N R f_0$")
small_ld_axis.set_ylabel("Small LD")
small_ld_axis.set_xlim([theory_xs[0],1e04])

# Plot DATA

n=1e05 # actual population size
fstars = numpy.logspace(-3,-1,20)
nstars = numpy.power(2,numpy.arange(3,11))
 
params = parameters.params
for type,symbol,counts_symbol in zip(['r'],['o'],['s']):
    
    sigmas = {fstar:[] for fstar in fstars}
    sigmasquareds = {fstar:[] for fstar in fstars}
    sigmafourths = {fstar:[] for fstar in fstars}
    les = {fstar:[] for fstar in fstars}
    small_lds = {fstar:[] for fstar in fstars}
    smaller_lds = {fstar:[] for fstar in fstars}
    denominatorsquareds = {fstar:[] for fstar in fstars}
    bare_numeratorsquareds = {fstar:[] for fstar in fstars}
    
    denominatorfourths = {fstar:[] for fstar in fstars}
    
    sigmasquareds_counts = {fstar:[] for fstar in fstars}
    
    gammas = []
    
    if debug:
        param_range = [10,20]
    else:
        param_range = range(0,len(params[type]))
    for param_idx in param_range:
    
        N = params[type][param_idx][2]
        r = params[type][param_idx][6]
        gamma = 2*N*r
        theta = 1.0/log(N)
        print "2Nr =", gamma
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

        n11s = numpy.array(n11s)*1.0
        n10s = numpy.array(n10s)*1.0
        n01s = numpy.array(n01s)*1.0
        ns = numpy.ones_like(n11s)*n
        n00s = ns-n10s-n11s-n01s
        
        smaller_lds = []
        smallerer_lds = []
        for m in nstars:
        
            fstar = 1.0/m
    
            print m
            good_200s = (n11s>1.5)*(n00s>(m-2.5))
            good_022s = (n10s>1.5)*(n01s>1.5)*(n00s>(m-4.5))
            good_111s = (n10s>0.5)*(n10s>0.5)*(n01s>0.5)*(n00s>(m-3.5))
    
            weights_200 = good_200s*n11s*(n11s-1)/2.0*numpy.exp(good_200s*(gammaln(n00s+1)-gammaln((n00s-(m-2))*good_200s+1)-gammaln(m-2+1)-gammaln(ns+1)+gammaln(ns-m+1)+gammaln(m+1)))
    
            print n11s[0:5]
            print n00s[0:5]
            print n10s[0:5]
            print n01s[0:5]
            print good_200s[0:5]

            print weights_200[0:5]
            
            weights_022 = good_022s*n10s*(n10s-1)/2.0*n01s*(n01s-1)/2.0*numpy.exp(good_022s*(gammaln(n00s+1)-gammaln((n00s-(m-4))*good_022s+1)-gammaln(m-4+1)-gammaln(ns+1)+gammaln(ns-m+1)+gammaln(m+1)))
    
            weights_111 = good_111s*n11s*n10s*n01s*numpy.exp(good_111s*(gammaln(n00s+1)-gammaln((n00s-(m-3))*good_111s+1)-gammaln(m-3+1)-gammaln(ns+1)+gammaln(ns-m+1)+gammaln(m+1)))
    
            P200 = weights_200.sum()
            P022 = weights_022.sum()
            P111 = weights_111.sum()
    
    
            smaller_ld = (m*m/2.0*P200-m/2.0*P111+P022)/(P022)/m
            
            smaller_lds.append(smaller_ld)
            
            small_ld_numerator = (numpy.square(f11s*f00s-f10s*f01s)*numpy.exp(-f11s/fstar-f10s/fstar-f01s/fstar)).mean()
            small_ld_denominator = (fAs*(1-fAs)*fBs*(1-fBs)*numpy.exp(-f11s/fstar-f10s/fstar-f01s/fstar)).mean()
            small_ld = small_ld_numerator*1.0/small_ld_denominator
            smallerer_lds.append(small_ld*m)
            
        collapse_xs = gamma/nstars
        collapse_ys = smaller_lds
        small_ld_axis.semilogx(collapse_xs, collapse_ys,'^', color='k',markersize=3)
        collapse_ys = smallerer_lds
        small_ld_axis.semilogx(collapse_xs, collapse_ys,'v', color='k',markersize=3)
        
        
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
            denominatorfourths[fstar].append(sigmafourth_denominator/theta**2)
            LE_numerator = (f11s*f10s*f01s*f00s*Hs).mean()
            LE_denominator = sigmafourth_denominator
            
            LE = LE_numerator/LE_denominator
            les[fstar].append(LE)
            
            small_ld_numerator = (f11s*f11s*numpy.exp(-f11s/fstar-f10s/fstar-f01s/fstar)).mean()
            small_ld_denominator = (f10s*f10s*f01s*f01s*numpy.exp(-f11s/fstar-f10s/fstar-f01s/fstar)).mean()
            small_ld = fstar*fstar*small_ld_numerator*1.0/small_ld_denominator
            small_lds[fstar].append(small_ld)
            
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
        les[fstar] = numpy.array(les[fstar])
        denominatorsquareds[fstar] = numpy.array(denominatorsquareds[fstar])
        bare_numeratorsquareds[fstar] = numpy.array(bare_numeratorsquareds[fstar])
        
        sigmasquareds_counts[fstar] = numpy.array(sigmasquareds_counts[fstar])
        
        denominatorfourths[fstar] = numpy.array(denominatorfourths[fstar])
        
        small_lds[fstar] = numpy.array(small_lds[fstar])
        
    for fstar in fstars:

        collapse_xs = gammas*fstar
        collapse_ys = sigmasquareds[fstar]/fstar
        
        colorVal = scalarMap.to_rgba(log10(fstar))
        
        
        line, = k2_axis.loglog(collapse_xs,collapse_ys,symbol,markersize=3,color=colorVal)
        
        
        
        #collapse_ys = sigmasquareds_counts[fstar]/fstar
        #pylab.loglog(collapse_xs,collapse_ys,counts_symbol,markersize=3)
        
        collapse_ys = sigmas[fstar]
        k1_axis.semilogx(collapse_xs, collapse_ys,symbol, color=colorVal,markersize=3)
        
        collapse_ys = les[fstar]
        le_axis.semilogx(collapse_xs, collapse_ys,symbol, color=colorVal,markersize=3)
        
        capped_fstar = min([fstar,1])
        collapse_ys = denominatorfourths[fstar]/capped_fstar/capped_fstar/capped_fstar/capped_fstar
        denom_axis.semilogx(collapse_xs, collapse_ys,symbol, color=colorVal,markersize=3)
        
        collapse_ys = small_lds[fstar]/fstar
        small_ld_axis.semilogx(collapse_xs, collapse_ys,symbol, color=colorVal,markersize=3)
        
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
k2_axis.semilogx(theory_xs,theory_ys,'k:')

    
theory_ys = numpy.array([ld_theory.scaled_neutral_ld(x) for x in theory_xs])
k2_axis.loglog(theory_xs, theory_ys,'k-',linewidth=1)

theory_ys = numpy.array([ld_theory.scaled_neutral_ld(x,k=1) for x in theory_xs])
k1_axis.semilogx(theory_xs, theory_ys,'k-',linewidth=1)

theory_ys = theory_xs
le_axis.semilogx(theory_xs, theory_ys,'k:')
le_axis.semilogx(theory_xs, numpy.ones_like(theory_ys),'k:')
le_axis.semilogx(theory_xs, numpy.zeros_like(theory_ys),'k:')

theory_ys = 1-2/theory_xs
le_axis.semilogx(theory_xs, theory_ys,'k:')

theory_ys = numpy.array([ld_theory.scaled_neutral_small_ld(x) for x in theory_xs])
small_ld_axis.loglog(theory_xs, theory_ys,'k-',linewidth=1)


k2_axis.set_ylim([3e-03,5])
k2_axis.set_xlim([3e-03,3e02])
k1_axis.set_xlim([3e-03,3e02])
k1_axis.set_ylim([-1,3])
k1_axis.set_yticks([-1,0,1,2,3])

m = k2_axis.scatter([1e-05],[1],c=[-2], vmin=vmin, vmax=vmax, cmap=cmap, marker='^')

f.subplots_adjust(top=0.89, hspace=0.05)
cax = f.add_axes([0.25, 0.95, 0.62, 0.02])
cbar = f.colorbar(m,cax=cax,orientation='horizontal',ticks=[-3,-2,-1])
#cbar.set_ticklabels(['$30$','$300$','$3000$'])
cbar.ax.tick_params(labelsize=8) 
f.text(0.115,0.94,'$\log_{10} f_0$')
#cbar.set_label('$U_d/U_b$')

le_axis.set_ylim([-0.1,1.2])
pylab.figure(1)
pylab.savefig('neutral_collapse.pdf',bbox_inches='tight')
pylab.figure(2)
if debug:
    filename = 'debug_neutral_le.pdf'
else:
    filename = 'neutral_le.pdf'
pylab.savefig(filename,bbox_inches='tight')

pylab.figure(3)
denom_axis.set_ylim([0,1.1])
pylab.savefig('neutral_denomiantor.pdf',bbox_inches='tight')

pylab.figure(4)
pylab.savefig('neutral_small_ld.pdf',bbox_inches='tight')

