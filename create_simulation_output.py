import sys
import numpy
import pylab
import gzip
import parameters
from math import log10,log
import ld_theory
from numpy.random import multinomial

type='smallr'
params = parameters.params
print "n11 n10 n01 n00 dl"
    
for param_idx in xrange(0,len(params[type])):

    rho = 2e-02
    
    N = params[type][param_idx][2]
    r = params[type][param_idx][6]
    ell = long(2*N*r/rho)
    n = N
    
    if ell>3e03:
        continue
    if ell<0.5:
        continue
        
    theta = 1.0/log(N)
    
    filename = 'output/output_%s_%d.txt.gz' % (type,param_idx)
    file = gzip.GzipFile(filename,"r")
    f11s = []
    f10s = []
    f01s = []
        
    # for counts (regenerated on the fly)
    n11s = []
    n10s = []
    n01s = []
    for line in file:
    
        if line.startswith('//'):
            continue
        items = line.split()
            
        f11 = float(items[0])
        f10 = float(items[1])
        f01 = float(items[2])
            
        n11 = long(f11*n)
        n10 = long(f10*n)
        n01 = long(f01*n)
        n00 = n-n11-n10-n01
    
        
        print "%d %d %d %d %d" % (n11,n10,n01,n00,ell)
                
    file.close()
   
