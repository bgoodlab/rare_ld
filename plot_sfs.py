import sys
import numpy
import pylab
import gzip
import parameters
from math import log10
import ld_theory

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
param_idx = 0
#for param_idx in xrange(0,len(params)):
if True:   
    N = params[param_idx][2]
    r = params[param_idx][6]
    Nr = N*r
       
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

print (2*fAs*(1-fAs)).mean()

bins = numpy.linspace(0,1,10000)
hist,edges = numpy.histogram(fAs,bins=bins)
hist = hist*1.0/hist.sum()
xs = bins[1:]
pylab.plot(xs,hist*xs,'r-')
hist,edges = numpy.histogram(fBs,bins=bins)
hist = hist*1.0/hist.sum()
pylab.semilogx(xs,hist*xs,'b-')
#pylab.xlim([0,0.01])
pylab.show()