import sys
import numpy
import pylab
import gzip
import parameters
from math import log10
import ld_theory

params = [(1e08,10,1e05,0,0,0,4e-04,1e-02)]

for param_idx in xrange(0,len(params)):
        
    filename = 'output_fstar_1.txt.gz'
    file = gzip.GzipFile(filename,"r")
    
    desired_num_samples = set([100000,300000,1000000,3000000,10000000,30000000,100000000])
    max_num_samples = max(desired_num_samples)
    sigmasquareds = []
    sigmafourths = []
    
    num_samples = 0
    running_numerator = 0
    running_numerator_squared = 0
    running_denominator = 0
    running_denominator_squared = 0
    print "Loading %s" % filename
    file.readline() # header
    for line in file:
    
        items = line.split()
        f11 = float(items[0])
        f10 = float(items[1])
        f01 = float(items[2])
        f00 = 1-f11-f10-f01
        
        fA = f11+f10
        fB = f11+f01
        
        D = f11*f00-f10*f01
        
        running_numerator+=D*D
        running_denominator+= (fA*(1-fA)*fB*(1-fB))
        
        running_numerator_squared += D**4
        running_denominator_squared+= (fA*(1-fA)*fB*(1-fB))**2
        
        num_samples+=1
        
        if num_samples in desired_num_samples:
            
            sigmasquared = running_numerator/running_denominator
            sigmasquareds.append(sigmasquared)
            
            sigmafourth = running_numerator_squared/running_denominator_squared/3.0
            sigmafourths.append(sigmafourth)
            
        if num_samples == max_num_samples:
            break
        
    file.close()
            
    sigmasquareds = numpy.array(sigmasquareds)
    sigmafourths = numpy.array(sigmafourths)
    
    desired_num_samples = numpy.array(sorted(desired_num_samples))
    
    etas = sigmafourths/sigmasquareds
    
    pylab.semilogx(desired_num_samples, etas,'k.-')
    theory_eta = ld_theory.neutral_eta(0.8)
    pylab.semilogx(desired_num_samples, numpy.ones_like(desired_num_samples)*theory_eta,'r:')
    pylab.savefig('fstar_test.pdf',bbox_inches='tight')
    