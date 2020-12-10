import sys
import gzip
import numpy
from numpy.random import multinomial
import parameters

params = parameters.params
n=1000

if __name__=='__main__':

    for type in ['eps','r']:
    
    
        for param_idx in [0]: #xrange(0,len(params[type])):
    
            input_filename = 'output/output_%s_%d.txt.gz' % (type,param_idx)
            file = gzip.GzipFile(input_filename,"r")
            
            output_filename = 'output/count_output_%s_%d.txt.gz' % (type,param_idx)
            output_file = gzip.GzipFile(output_filename,"w")
            
            print "Calculating %s" % output_filename
            for line in file:
    
                if line.startswith('//'):
                    output_file.write(line)
                    continue
                items = line.split()
            
                f11 = float(items[0])
                f10 = float(items[1])
                f01 = float(items[2])
                f00 = 1-f11-f10-f01
                f00 = f00*(f00>0)
            
                fs = numpy.array([f11,f10,f01,f00])
                fs = fs/fs.sum()
            
            
                ns = multinomial(n,fs)
            
                n11 = ns[0]
                n10 = ns[1]
                n01 = ns[2]
                n00 = n-n11-n10-n01
                
                output_file.write("%g %g %g %d %d %d %d" % (f11,f10,f01,n11,n10,n01,n00))
                
            file.close()
            output_file.close()