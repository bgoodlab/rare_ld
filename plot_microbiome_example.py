import sys
import numpy
import pylab
import gzip
from scipy.special import gammaln 

import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
import matplotlib.gridspec as gridspec
from numpy.random import randint, shuffle, poisson, binomial, choice, hypergeometric
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

def calculate_sigmasquared(n11s,n10s,n01s,n00s,ntots,fstar):
    
    # Now do version based on counts
    # First calculate numerator
    rsquared_numerators = n11s*(n11s-1)*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-0)          
    rsquared_numerators += -2*n10s*n01s*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)
    rsquared_numerators += n10s*(n10s-1)*n01s*(n01s-1)*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-2)
    
    # Divide by sample size
    rsquared_numerators = rsquared_numerators*(ntots>3.5)*1.0/(ntots*(ntots-1)*(ntots-2)*(ntots-3)+10*(ntots<3.5))
    
    #1
    rsquared_denominators = n10s*(n10s-1)*n01s*(n01s-1)*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-2)
    #2
    rsquared_denominators += n10s*n01s*(n01s-1)*n00s*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-2)       
    #3
    rsquared_denominators += n10s*(n10s-1)*n01s*n11s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-1)
    #4
    rsquared_denominators += n10s*n01s*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)           
    #5
    rsquared_denominators += n10s*(n10s-1)*n01s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-1)       
    #6
    rsquared_denominators += n10s*n01s*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)
    #7
    rsquared_denominators += n10s*(n10s-1)*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-0)
    #8
    rsquared_denominators += n10s*n11s*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-0)
    #9
    rsquared_denominators += n10s*n01s*(n01s-1)*n11s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-2)
    #10
    rsquared_denominators += n01s*(n01s-1)*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-2)
    #11
    rsquared_denominators += n10s*n01s*n11s*(n11s-1)*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)
    #12
    rsquared_denominators += n01s*n11s*(n11s-1)*n00s*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-1)
    #13
    rsquared_denominators += n10s*n01s*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)
    #14
    rsquared_denominators += n01s*n11s*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-1)
    #15
    rsquared_denominators += n10s*n11s*(n11s-1)*n00s*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-0)
    #16
    rsquared_denominators += n11s*(n11s-1)*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-0)
    
    # divide by sample size
    rsquared_denominators = rsquared_denominators*(ntots>3.5)*1.0/(ntots*(ntots-1)*(ntots-2)*(ntots-3)+10*(ntots<3.5))

    return rsquared_numerators.mean()/rsquared_denominators.mean()
    
debug = True

species_name = sys.argv[1].split("_")[0]

mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

theory_xs = numpy.logspace(-3,3.5,100)

# Set up figure
f, axes = pylab.subplots(2,2,squeeze=False)
f.set_size_inches(6,4.5)

ldl_axis = axes[0][0]
ldl_axis.set_ylabel("LD, $\\sigma_d^2(f_0)$")
ldl_axis.set_xlabel("Distance between SNVs, $\ell$")
ldl_axis.set_ylim([6e-04,6e-01])
ldl_axis.set_xlim([2,8e03])

ldl_axis.fill_between([4e03,1e04],[1,1],[1e-04,1e-04],color='0.8')

ldl_axis.loglog([1e-06],[1e-06],'k.')

ldf_axis = axes[0][1]
#ldf_axis.set_ylabel("$\\sigma_d^2")
ldf_axis.set_xlabel("Frequency scale, $f_0$")
ldf_axis.set_ylim([6e-04,6e-01])
ldf_axis.set_xlim([2e-02,2])
ldf_axis.loglog([1e-06],[1e-06],'k.')
ldf_axis.fill_between([1.0,2.0],[1,1],[1e-04,1e-04],color='0.8')

sfs_axis = axes[1][0]
sfs_axis.set_ylabel("Fraction SNV pairs")
sfs_axis.set_xlabel("Single mutant count, $n_A$")

n11_axis = axes[1][1]
n11_axis.set_xlabel("Double mutant count, $n_{AB}$")

pylab.subplots_adjust(hspace=0.35,wspace=0.25)


#fstars = [1e02,1e-01]
fstars = [1e02,1e-01,3e-02]
fstar_colors = ['#045a8d','#2b8cbe','#74a9cf']
fstar_labels=['$f_0=\\infty$','$f_0=0.1$', '$f_0=0.03$']
#file = gzip.GzipFile("fragilis_ben_ld.txt.gz","r")
#file = gzip.GzipFile("shahii_ben2_ld.txt.gz","r")
#file = gzip.GzipFile("rectale_ben_ld.txt.gz","r")
#file = gzip.GzipFile("stercoris_ben_ld.txt.gz","r")
file = gzip.GzipFile(sys.argv[1],"r")


file.readline() # header

n11s = []
n10s = []
n01s = []
n00s = []
ells = []
print "Loading data..."
for line in file:
    
    items = line.split()
    n11 = long(items[0])
    n10 = long(items[1])
    n01 = long(items[2])
    n00 = long(items[3])
    ell = long(items[4])

    n11s.append(n11)
    n10s.append(n10)
    n01s.append(n01)
    n00s.append(n00)
    ells.append(ell)

file.close()

n11s =numpy.array(n11s)*1.0
n10s = numpy.array(n10s)*1.0
n01s = numpy.array(n01s)*1.0
n00s = numpy.array(n00s)*1.0
ntots = n11s+n10s+n01s+n00s
ells = numpy.array(ells)

max_ntot = ntots.max()

good_idxs = (ntots>0.95*max_ntot)
print "max_ntot", max_ntot
print "Kept", good_idxs.sum()*1.0/len(good_idxs)

n11s = n11s[good_idxs]
n10s = n10s[good_idxs]
n01s = n01s[good_idxs]
n00s = n00s[good_idxs]
ntots = ntots[good_idxs]
ells = ells[good_idxs]

min_ntot = ntots.min()

f11s = n11s*1.0/ntots
f10s = n10s*1.0/ntots
f01s = n01s*1.0/ntots
f00s = n00s*1.0/ntots 

fAs = f11s+f10s
fBs = f11s+f01s   

nAs = n11s+n10s
nBs = n11s+n10s


macs = numpy.fmin(nAs,ntots-nAs)
bins = numpy.arange(0,max_ntot)+0.5
h,bin_edges = numpy.histogram(macs,bins=bins)
ks = numpy.arange(1,max_ntot)

minor_h = (h+h[::-1])*1.0
minor_h /= minor_h.sum()
minor_h *= 2

fmin = 0.13
fmax = 0.17


sfs_axis.loglog(ks,minor_h,label='Observed')
sfs_axis.set_xlim([1,min_ntot/2.0])

nstar = long(min_ntot/2.0)
ylims = sfs_axis.get_ylim()
bottom = ylims[0]
top = ylims[1]

print bottom, top

#sfs_axis.fill_between([fmin*min_ntot,fmax*min_ntot],[top,top],[bottom,bottom],color='0.7')

print h[0:5]
print ks[0:5]
#sfs_axis.loglog([2,4],[h[0],h[0]/(5.0/2)],'k:')
#sfs_axis.loglog([1.3,4],[h[0]/3,h[0]/3/(4.0/1.3)**2],'k:')

neutral_h = 1.0/(ks*(max_ntot-ks))


neutral_h = neutral_h*minor_h[nstar-1]/neutral_h[nstar-1]
sfs_axis.loglog(ks,neutral_h,'k:',label='$1/n_A(n-n_A)$')

sfs_axis.set_ylim([bottom,top])

sfs_axis.legend(loc='lower left',frameon=False,numpoints=1)
ellranges = [(1e04,1e08),(1e03,2e03),(1e02,3e02)]
colors = ['#1f77b4','#2ca02c','#ff7f0e']
labels =['Genome avg','1000<$\ell$<2000','100<$\ell$<300']

for idx in xrange(0,len(ellranges)):
    ellmin,ellmax = ellranges[idx]
    color = colors[idx]
    
    bins = numpy.arange(0,ntots.max())-0.5
    xs = numpy.arange(0,ntots.max()-1)
    
    good_idxs = (ells>ellmin)*(ells<ellmax)*(fAs<fmax)*(fBs<fmax)*(fAs>fmin)*(fBs>fmin)

    h,bin_edges = numpy.histogram(n11s[good_idxs],bins=bins)
    h = h*1.0/h.sum()
    n11_axis.semilogy(xs,h,'.-',markersize=3,color=color)

n11_axis.set_xlim([-0.5,20])

ylims = n11_axis.get_ylim()
bottom = ylims[0]
top = ylims[1]
n11_axis.plot([min_ntot*fmin*fmin,max_ntot*fmin*fmin],[bottom,top],'k:')
n11_axis.plot([min_ntot*fmin,min_ntot*fmin],[bottom,top],'k:')

n11_axis.set_ylim([bottom,top])


f0s = numpy.hstack([numpy.logspace(-1.5,-0.5,20),[1e02]])
capped_f0s = numpy.clip(f0s,0,1.5)


for idx in xrange(0,len(ellranges)):
    ellmin,ellmax = ellranges[idx]
    color = colors[idx]
    label=labels[idx]
    
    bins = numpy.arange(0,ntots.max())-0.5
    xs = numpy.arange(0,ntots.max()-1)
    
    good_idxs = (ells>ellmin)*(ells<ellmax)

    sigmasquareds = []
    for f0 in f0s:
        print f0
        sigmasquared = calculate_sigmasquared(n11s[good_idxs],n10s[good_idxs],n01s[good_idxs],n00s[good_idxs],ntots[good_idxs],f0)
        sigmasquareds.append(sigmasquared)
    
    sigmasquareds = numpy.array(sigmasquareds)    
    line, = ldf_axis.loglog(capped_f0s[:-1],sigmasquareds[:-1],'-',markersize=3,color=color,label=label)
    ldf_axis.loglog(capped_f0s[-2:],sigmasquareds[-2:],'.:',color=color,markersize=3)

ldf_axis.legend(loc='upper left',frameon=False,numpoints=1)
pylab.savefig('%s_example.pdf' % species_name,bbox_inches='tight')
 
raw_theory_ells = numpy.logspace(0,3,200)*3
theory_ells = numpy.hstack([raw_theory_ells,[6e03]])
big_theory_ells = numpy.hstack([raw_theory_ells,[1e07]])

print "Done!"
for idx in xrange(0,len(fstars)):
    fstar = fstars[idx]
    color=fstar_colors[idx]
    label=fstar_labels[idx]
    
    capped_fstar = min([fstar,1.0])
    print "Processing fstar =", fstar, capped_fstar
    
    # Now do version based on counts
    # First calculate numerator
    rsquared_numerators = n11s*(n11s-1)*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-0)          
    rsquared_numerators += -2*n10s*n01s*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)
    rsquared_numerators += n10s*(n10s-1)*n01s*(n01s-1)*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-2)
    
    # Divide by sample size
    rsquared_numerators = rsquared_numerators*(ntots>3.5)*1.0/(ntots*(ntots-1)*(ntots-2)*(ntots-3)+10*(ntots<3.5))
    
    #1
    rsquared_denominators = n10s*(n10s-1)*n01s*(n01s-1)*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-2)
    #2
    rsquared_denominators += n10s*n01s*(n01s-1)*n00s*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-2)       
    #3
    rsquared_denominators += n10s*(n10s-1)*n01s*n11s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-1)
    #4
    rsquared_denominators += n10s*n01s*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)           
    #5
    rsquared_denominators += n10s*(n10s-1)*n01s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-1)       
    #6
    rsquared_denominators += n10s*n01s*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-0)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)
    #7
    rsquared_denominators += n10s*(n10s-1)*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-2+n01s-0)
    #8
    rsquared_denominators += n10s*n11s*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-0)
    #9
    rsquared_denominators += n10s*n01s*(n01s-1)*n11s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-2)
    #10
    rsquared_denominators += n01s*(n01s-1)*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-2)
    #11
    rsquared_denominators += n10s*n01s*n11s*(n11s-1)*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)
    #12
    rsquared_denominators += n01s*n11s*(n11s-1)*n00s*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-1)
    #13
    rsquared_denominators += n10s*n01s*n11s*n00s*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-1)
    #14
    rsquared_denominators += n01s*n11s*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-1)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-1)
    #15
    rsquared_denominators += n10s*n11s*(n11s-1)*n00s*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-1+n01s-0)
    #16
    rsquared_denominators += n11s*(n11s-1)*n00s*(n00s-1)*numpy.power(1-2.0/ntots/fstar,n11s-2)*numpy.power(1-1.0/ntots/fstar,n10s-0+n01s-0)
    
    # divide by sample size
    rsquared_denominators = rsquared_denominators*(ntots>3.5)*1.0/(ntots*(ntots-1)*(ntots-2)*(ntots-3)+10*(ntots<3.5))

    sigmasquareds2 = []
    avg_ells2 = []
    for l in big_theory_ells:
        
        lmin = l*numpy.power(10,-0.1)
        lmax = l*numpy.power(10,0.1)
        
        good_idxs = (ells>=lmin)*(ells<=lmax)
        
        if good_idxs.sum() < 10:
            # fewer than 10 points!
            continue
        
        sigmasquared = (rsquared_numerators*(ells>=lmin)*(ells<=lmax)).mean() / (rsquared_denominators*(ells>=lmin)*(ells<=lmax)).mean()
       
        local_ells = ells[(ells>=lmin)*(ells<=lmax)]
        
        # Harmonic mean
        avg_ell = 1.0/(1.0/local_ells).mean()
        # regular mean
        # avg_ell = local_ells.mean()
        
        avg_ells2.append(avg_ell)
        sigmasquareds2.append(sigmasquared)
        
        #print l,lmin,lmax, avg_ell, sigmasquared
        
        
    
    avg_ells = numpy.array(avg_ells2)
    avg_ells[-1]=6e03
    sigmasquareds = numpy.array(sigmasquareds2)
    
    line, = ldl_axis.loglog(avg_ells[:-1],sigmasquareds[:-1],'-',label=label,color=color)
    ldl_axis.loglog(avg_ells[-2:],sigmasquareds[-2:],'.:',color=color,markersize=3)

ldl_axis.text(6e03,1.7e-04,'Genome\navg',     horizontalalignment='center',fontsize='7')
    
ldl_axis.legend(loc='lower left',frameon=False,numpoints=1)        
#theory_sigmasquareds = numpy.array([neutral_rsquared(ell*0.005) for ell in theory_ells])
#ldl_axis.loglog(theory_ells[:-1],theory_sigmasquareds[:-1],'k-',linewidth=0.5) 
#pylab.ylim([1e-02,1])

pylab.savefig('%s_example.pdf' % species_name,bbox_inches='tight')
