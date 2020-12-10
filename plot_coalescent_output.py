import pylab
import numpy
import gzip
import sys
from math import log
import ld_theory

filename = sys.argv[1]
file = gzip.GzipFile(filename,"r")

# Process data

fABs = []
fAs = []
fBs = []
double_weights = []
single_weights = []
Ttots = []

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
    
    fAB = nAB*1.0/n
    fAb = nAb*1.0/n
    faB = naB*1.0/n
    fab = nab*1.0/n
    
    fA = fAB+fAb
    fB = fAB+faB
    
    weights_x.append(Tx)
    weights_y.append(Ty)
    
    single_weights.append(Tx)
    double_weights.append(Tx*Ty)
    fABs.append(fAB)
    fAs.append(fA)
    fBs.append(fB)
    Ttots.append(Tx)
    Ttots.append(Ty)
    
fABs = numpy.array(fABs)
fAs = numpy.array(fAs)
fBs = numpy.array(fBs)
double_weights = numpy.array(double_weights)
single_weights = numpy.array(single_weights)
Ttots = numpy.array(Ttots)

weights_x = numpy.array(weights_x)
weights_y = numpy.array(weights_y)

print weights_x.mean()/log(1000)
print weights_y.mean()/log(1000)
pi = (2*fAs*(1-fAs)*single_weights).mean()
print "Pi=",pi
print "Pi_W=", single_weights.mean()/log(1000)
sigmasquareds = []
fstars = numpy.logspace(-2,0,20)
for fstar in fstars:
    
    numerator = (numpy.square(fABs-fAs*fBs)*double_weights*numpy.exp(-fAs/fstar-fBs/fstar)).mean()
    denominator = (fAs* (1-fAs)*fBs*(1-fBs)*double_weights*numpy.exp(-fAs/fstar-fBs/fstar)).mean()

    sigmasquared = numerator*1.0/denominator
    
    sigmasquareds.append(sigmasquared)

sigmasquareds = numpy.array(sigmasquareds)

Nr=10
rhos = 2*Nr*fstars

pylab.loglog(rhos,sigmasquareds/fstars,'b.')

theory_ys = numpy.array([ld_theory.scaled_neutral_ld(rho) for rho in rhos])
pylab.loglog(rhos, theory_ys,'r-')

theory_ys = numpy.array([ld_theory.scaled_asexual_sigmasquared(rho) for rho in rhos])
pylab.loglog(rhos, theory_ys,'g-')


#pylab.loglog(fstars,5.0/11*numpy.ones_like(fstars),'r:')
#pylab.loglog(fstars,fstars,'r:')
#pylab.loglog(fstars,fstars*numpy.log(1/fstars+1)*0.5,'r:')

pylab.savefig('coalescent_fig.pdf',bbox_inches='tight')