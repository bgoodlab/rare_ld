import numpy
import sys

params = {}

num_runs = 100000000
N = 1e05
dt = 100
s1 = 0
s2 = 0
eps = 0
rs = numpy.hstack([[0],numpy.logspace(-6,-1,21)])

# type == 'r'
params['r'] = [(num_runs,dt,N,s1,s2,eps,r) for r in rs]

ss = numpy.logspace(-5,-2,13)
params['selA'] = [(num_runs,dt,N,s,0,0,0) for s in ss]
params['selB'] = [(num_runs,dt,N,0,s,0,0) for s in ss]
params['selAB'] = [(num_runs,dt,N,s/2,s/2,0,0) for s in ss]
params['eps'] = [(num_runs,dt,N,0,0,s,0) for s in ss]
params['negeps'] = [(num_runs,dt,N,s,s,-s,0) for s in ss]

s1 = 1e-03/8
s2 = 1e-03/8
eps = 0 
rs = numpy.hstack([[0],numpy.logspace(-6,-2,17)])
params['r_selA'] = [(num_runs,dt,N,s1,s2,eps,r) for r in rs]

N = 3e03
num_runs = 1000000
dt = 10
s1 = 0
s2 = 0
eps = 0
rs = numpy.hstack([[0],numpy.logspace(-5,-1,15)])
params['smallr'] = [(num_runs,dt,N,s1,s2,eps,r) for r in rs]


if __name__=='__main__':

    if sys.argv[1]=='idxs':
        type = sys.argv[2].strip()
        for r_idx in xrange(0,len(params[type])):
            print r_idx
    elif sys.argv[1]=='get_params':
        idx = long(sys.argv[3])
        type = sys.argv[2].strip()
        print " ".join([str(item) for item in params[type][idx]])