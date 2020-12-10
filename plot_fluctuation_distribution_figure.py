import sys
import numpy
import pylab
import gzip
import parameters
from math import log10,log
import ld_theory
from numpy.random import multinomial,binomial

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

theory_xs = numpy.logspace(-3,3.5,50)

# Set up figure
pylab.figure(1,figsize=(3.42,2.5))
f = pylab.gcf()

n11_axis = pylab.gca()
n11_axis.set_ylabel("Probability, $p(f_{AB}|f_{A}, f_{B} \\approx f_0)$")
n11_axis.set_xlabel("Double mutant frequency, $f_{AB}$")

focal_idxs = [5,9,12]
focal_colors = ['#ff7f0e','#2ca02c','#1f77b4']
focal_labels = ['$\\rho \\ll 1/f_0$','$\\rho < 1/f_0$','$\\rho > 1/f_0$']
fmin=0.13
fmax=0.17
small_n=2e02

n11_axis.set_xlim([-0.01,0.2])
n11_axis.set_ylim([1e-03,1])

n11_axis.semilogy([fmin*fmin,fmin*fmin],[1e-04,1],'k:')
n11_axis.semilogy([fmin,fmin],[1e-04,0.1],'k:')


#n11_axis.fill_between([fmin*fmin,fmax*fmax],[1,1],[1e-04,1e-04],color='0.8')
#n11_axis.fill_between([fmin,fmax],[1,1],[1e-04,1e-04],color='0.8')

# Plot DATA

n=1e05 # actual population size
small_n=2e02
#fstars = numpy.logspace(-3,-1,20)
fstars = numpy.array([0.001,0.01,0.03,0.1])
params = parameters.params
for type,symbol,counts_symbol in zip(['r'],['o'],['s']):
    
    for param_idx,color,other_label in zip(focal_idxs,focal_colors,focal_labels):
    
        N = params[type][param_idx][2]
        r = params[type][param_idx][6]
        gamma = 2*N*r
        theta = 1.0/log(N)
    
        #label=('$\\rho=%0.1f$, $\\rho f_0=%0.1f$' % (gamma*fmin,gamma*fmin*fmin))
        
        label=('$\\rho=%0.1f$  (%s)' % (gamma*fmin,other_label))
    
        
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
            
            fA = f11+f10
            fB = f11+f01
            
            if not (fA>fmin)*(fA<fmax)*(fB>fmin)*(fB<fmax):
                continue
            
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
        
        #small_n11s = binomial(small_n,f11s)
        small_n11s = numpy.around(f11s*small_n)
        
        bins = numpy.arange(0,small_n)-0.5
        xs = numpy.arange(0,small_n-1)*1.0/small_n
    
        h,bin_edges = numpy.histogram(small_n11s,bins=bins)
        n11_axis.semilogy(xs,h*1.0/h.sum(),'.-',markersize=3,color=color,label=label)

        print "Loaded", len(ns), "pairs!"
    
    
    

n11_axis.legend(frameon=False,loc='upper center',numpoints=1,scatterpoints=1)

pylab.savefig('neutral_fluctuation_distribution.pdf',bbox_inches='tight')
