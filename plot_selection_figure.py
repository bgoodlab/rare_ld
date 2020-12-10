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
k2_axis.set_xlabel("$2 N s f_0$")
k2_axis.set_xlim([theory_xs[0],theory_xs[-1]])

vmin=-3
vmax=-1
cmap='jet_r'

jet = cm = pylab.get_cmap(cmap) 
cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)



# Plot DATA

n=1e05 # actual population size
fstars = numpy.logspace(-3,-1,10)
#fstars = numpy.array([1e-03,1e-02,1e-01])
params = parameters.params
for type,symbol,counts_symbol in zip(['selAB','negeps','eps'],['o','^','v'],['s','s','s']):
    
    sigmas = {fstar:[] for fstar in fstars}
    sigmasquareds = {fstar:[] for fstar in fstars}
    sigmafourths = {fstar:[] for fstar in fstars}
    denominatorsquareds = {fstar:[] for fstar in fstars}
    bare_numeratorsquareds = {fstar:[] for fstar in fstars}
    
    sigmasquareds_counts = {fstar:[] for fstar in fstars}
    
    gammas = []
    
    for param_idx in xrange(0,len(params[type])):
    
        N = params[type][param_idx][2]
        sA = params[type][param_idx][3]
        sB = params[type][param_idx][4]
        eps = params[type][param_idx][5]
        gamma = 2*(N*eps+N*sA+N*sB)
        theta = 1.0/log(N)
    
        thetaA = ld_theory.calculate_theta(N,sA)
        thetaB = ld_theory.calculate_theta(N,sB)
        print N,sA,sB,eps,gamma
    
        if gamma<20:
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
            
            denominatorsquareds[fstar].append(sigmasquared_denominator/thetaA/thetaB)
            
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
            
            bare_numeratorsquareds[fstar].append(( numpy.square(bare_Ds)*Hs).mean()/thetaA/thetaB)

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
        
        if type=='eps':
            good_idxs = (collapse_xs<1.0/numpy.sqrt(fstar))
        else:
            good_idxs = (collapse_xs<1e08)    
            
        
        colorVal = scalarMap.to_rgba(log10(fstar))
        
        line, = k2_axis.loglog(collapse_xs[good_idxs],collapse_ys[good_idxs],symbol,markersize=3,color=colorVal)
        
        #collapse_ys = sigmasquareds_counts[fstar]/fstar
        #pylab.loglog(collapse_xs,collapse_ys,counts_symbol,markersize=3)
        
        collapse_ys = sigmas[fstar]
        k1_axis.semilogx(collapse_xs[good_idxs], collapse_ys[good_idxs],symbol, color=colorVal,markersize=3)
        
        
        #pylab.figure(2)
        #collapse_ys = sigmafourths[fstar]/sigmasquareds[fstar]
        #pylab.loglog(collapse_xs,collapse_ys,symbol,color=color,markersize=3)
        
        #pylab.figure(3)
        #collapse_ys = denominatorsquareds[fstar]/fstar**2
        #pylab.semilogx(collapse_xs,collapse_ys,symbol,color=color,markersize=3)
        
        #pylab.figure(4)
        #collapse_ys = bare_numeratorsquareds[fstar]/fstar**3
        #pylab.loglog(collapse_xs,collapse_ys,symbol,color=color,markersize=3)
        
        
    
theory_ys = numpy.array([ld_theory.scaled_asexual_ld(x,0,0,k=2) for x in theory_xs])
#k2_axis.loglog(theory_xs, theory_ys,'k-',linewidth=1)
theory_ys = numpy.array([ld_theory.scaled_asexual_ld(x/2,x/2,0,k=2) for x in theory_xs])
k2_axis.loglog(theory_xs, theory_ys,'k-',linewidth=1)
theory_ys = numpy.array([ld_theory.scaled_asexual_ld(0,0,x,k=2) for x in theory_xs])
k2_axis.loglog(theory_xs, theory_ys,'k-',linewidth=1)
theory_ys = numpy.array([ld_theory.scaled_asexual_ld(x,x,-x,k=2) for x in theory_xs])
k2_axis.semilogx(theory_xs, theory_ys,'k-',linewidth=1)

theory_ys = numpy.array([ld_theory.scaled_asexual_ld(x,0,0,k=1) for x in theory_xs])
#k1_axis.semilogx(theory_xs, theory_ys,'k-',linewidth=1)
theory_ys = numpy.array([ld_theory.scaled_asexual_ld(x/2,x/2,0,k=1) for x in theory_xs])
k1_axis.semilogx(theory_xs, theory_ys,'k-',linewidth=1)
theory_ys = numpy.array([ld_theory.scaled_asexual_ld(0,0,x,k=1) for x in theory_xs])
k1_axis.semilogx(theory_xs, theory_ys,'k-',linewidth=1)
theory_ys = numpy.array([ld_theory.scaled_asexual_ld(x,x,-x,k=1) for x in theory_xs])
k1_axis.semilogx(theory_xs, theory_ys,'k-',linewidth=1)

theory_ys = 2.0/numpy.square(theory_xs)
#k2_axis.plot(theory_xs,theory_ys,'k:')
k2_axis.set_ylim([3e-03,5])
k2_axis.set_xlim([3e-03,3e02])
k1_axis.set_xlim([3e-03,3e02])
k1_axis.set_ylim([-1,4])
k1_axis.set_yticks([-1,0,1,2,3])

m = k2_axis.scatter([1e-05],[1],c=[-2], vmin=vmin, vmax=vmax, cmap=cmap, marker='^')

k2_axis.plot([1e-05],[1e-05],'ko',label='$s_A=s_B=s/2, \, \epsilon=0$',markersize=3)
k2_axis.plot([1e-05],[1e-05],'k^',label='$s_A=s_B=s, \, \epsilon=-s$',markersize=3)
k2_axis.plot([1e-05],[1e-05],'kv',label='$s_A=s_B=0, \, \epsilon=+s$',markersize=3)
#k2_axis.plot([1e-05],[1e-05],'k<',label='$s_A=s, \, s_B=\epsilon=0$',markersize=3)

k2_axis.legend(frameon=False,loc='lower left',numpoints=1,scatterpoints=1)

f.subplots_adjust(top=0.89, hspace=0.05)
cax = f.add_axes([0.25, 0.95, 0.62, 0.02])
cbar = f.colorbar(m,cax=cax,orientation='horizontal',ticks=[-3,-2,-1])
#cbar.set_ticklabels(['$30$','$300$','$3000$'])
cbar.ax.tick_params(labelsize=8) 
f.text(0.115,0.94,'$\log_{10} f_0$')
#cbar.set_label('$U_d/U_b$')


pylab.savefig('selection_collapse.pdf',bbox_inches='tight')
