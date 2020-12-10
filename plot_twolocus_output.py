import sys
import numpy
import pylab
import gzip

fstars = numpy.logspace(-3,0,10)
pylab.loglog(fstars, numpy.ones_like(fstars)*5.0/11,'k:')
pylab.loglog(fstars, 5.0/11*numpy.sqrt(fstars),'k:')
pylab.loglog(fstars, 5.0/11*fstars,'k:')
pylab.loglog(fstars, 0.5*fstars*numpy.log(1.0/fstars+1),'r:')


Nss = [0,10,100]
for Ns in Nss:
    filename = "output_%d.txt.gz" % Ns
    file = gzip.GzipFile(filename,"r")
    f11s = []
    f10s = []
    f01s = []

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

    Ds = f11s*f00s-f10s*f01s
    Dsquareds = numpy.square(Ds)
    sampling_variances = fAs*(1-fAs)*fBs*(1-fBs)


    sigmasquared_numerators = []
    sigmasquared_denominators = []

    for fstar in fstars:
    
        good_idxs = (fAs<=fstar)*(fBs<=fstar)
    
        #sigmasquared_numerators.append( Dsquareds[good_idxs].mean() )
        #sigmasquared_denominators.append( sampling_variances[good_idxs].mean() )

        sigmasquared_numerators.append( (Dsquareds*numpy.exp(-fAs/fstar-fBs/fstar)).mean() )
        sigmasquared_denominators.append( (sampling_variances*numpy.exp(-fAs/fstar-fBs/fstar)).mean() )


    sigmasquared_numerators = numpy.array(sigmasquared_numerators)
    sigmasquared_denominators = numpy.array(sigmasquared_denominators)
    sigmasquareds = sigmasquared_numerators/sigmasquared_denominators


    print sigmasquareds[-1]
    pylab.loglog(fstars, sigmasquareds, '.-')


#frequency_bins = numpy.linspace(0,1,1000)        
#pAs = numpy.histogram(fAs,bins=frequency_bins)[0]
#pBs = numpy.histogram(fBs,bins=frequency_bins)[0]
   
#print len(pAs), len(frequency_bins)  
#pylab.loglog(frequency_bins[1:], pAs,'-')
#pylab.loglog(frequency_bins[1:], pBs,'-')

pylab.savefig('twolocus_output.pdf',bbox_inches='tight')