import numpy
from numpy.random import choice, randint, exponential,uniform
from scipy.misc import comb
import sys
from math import log

num_tree_replicates = 1000
num_mutation_replicates = 100

sample_size = 1000

Nr = float(sys.argv[1])

blank_chromosome = frozenset()
ancestral_chromosome = frozenset(numpy.arange(0,sample_size))

avg_TtotA = 0
avg_TtotB = 0

for tree_idx in xrange(0,num_tree_replicates):

    # set up sample
    sample = [(frozenset((i,)),frozenset((i,))) for i in xrange(0,sample_size)]
    
    lineage_times_a = {lineage[0]:0 for lineage in sample}
    lineage_times_b = {lineage[1]:0 for lineage in sample}
    
    # Adding these entries will make main loop simpler
    # we will remove them at end
    lineage_times_a[blank_chromosome]=0
    lineage_times_b[blank_chromosome]=0
    
    n = len(sample)
    while True:
        
        # Calculate recombination eligible ones
        # (save time by not recombining things you don't need!)
        recombination_eligible=[]
        if Nr>0:
            for i in xrange(0,len(sample)):
            
                if len(sample[i][0])==0 or len(sample[i][1])==0:
                    continue                
                else:
                    recombination_eligible.append(i)
        nr = len(recombination_eligible)
        n = len(sample)   
        # Sample time to next event
        coalescent_rate = n*(n-1)/2.0
        recombination_rate = Nr*nr
        total_rate = coalescent_rate+recombination_rate
        
        recombination_probability = recombination_rate/(total_rate)
        
        T = exponential(1.0/total_rate)
        
        #print sample
        #print recombination_rate
        #print coalescent_rate
        #print 1.0/total_rate, T
        
        #print "Recording", n, T
        for lineage in sample:
            lineage_times_a[lineage[0]]+=T
            lineage_times_b[lineage[1]]+=T
        
        if uniform() < recombination_probability:
            # A recombination event!
            
            # choose a random individual to recombine
            i = recombination_eligible[randint(0,nr)] 
                
            new_lineage_1 = (sample[i][0],blank_chromosome)
            new_lineage_2 = (blank_chromosome,sample[i][1])
                
            sample[i] = new_lineage_1
            sample.append(new_lineage_2)
                
            # both lineages are already present on their respective chromosomes
            # so don't need to add anything to lineage times
            
        else:
            if n==2:
                break
                
            # A coalescent event!
            # Sample individuals to coalesce
            # (random pair)
            i = randint(0,n)
            j = randint(0,n-1)
            if j>=i:
                j+=1
        
            #print "Merging", i,j
            # merge j with i
            new_lineage = ((sample[i][0]|sample[j][0]),(sample[i][1]|sample[j][1]))
            
            if new_lineage[0] not in lineage_times_a:
                lineage_times_a[new_lineage[0]]=0
            
            if new_lineage[1] not in lineage_times_b:
                lineage_times_b[new_lineage[1]]=0
            
            sample[i] = new_lineage
            del sample[j]
        
        
        n = len(sample)
        #print n
        
    # Done sampling tree
    
    # Delete blank chromosomes
    # (not in present day sample)
    del lineage_times_a[blank_chromosome]
    if ancestral_chromosome in lineage_times_a:
        del lineage_times_a[ancestral_chromosome]
    del lineage_times_b[blank_chromosome]
    if ancestral_chromosome in lineage_times_b:
        del lineage_times_b[ancestral_chromosome]
    
    #print lineage_times_a
    
    # Calculate times for each lineage
    lineages_a = lineage_times_a.keys()
    lineage_idxs_a = numpy.arange(0,len(lineages_a))
    Ts = numpy.array([lineage_times_a[lineage] for lineage in lineages_a])
    Ttot_a = Ts.sum()*1.0
    ps_a = Ts/Ttot_a
    
    lineages_b = lineage_times_b.keys()
    lineage_idxs_b = numpy.arange(0,len(lineages_b))
    Ts = numpy.array([lineage_times_b[lineage] for lineage in lineages_b])
    Ttot_b = Ts.sum()*1.0
    ps_b = Ts/Ttot_b
    
    #print Ttot
    mutation_branches_a = choice(lineage_idxs_a,size=num_mutation_replicates,p=ps_a)
    mutation_branches_b = choice(lineage_idxs_b,size=num_mutation_replicates,p=ps_b)
    
    avg_TtotA += Ttot_a
    avg_TtotB += Ttot_b
    
    #print mutation_branches[0]
    for mutation_idx in xrange(0,num_mutation_replicates):
        nA = len(lineages_a[mutation_branches_a[mutation_idx]])
        nB = len(lineages_b[mutation_branches_b[mutation_idx]])
        nAB = len(lineages_a[mutation_branches_a[mutation_idx]] & lineages_b[mutation_branches_b[mutation_idx]])
        
        nAb = nA-nAB
        naB = nB-nAB
        nab = sample_size-nAB-nAb-naB
        
        print nAB,nAb,naB,nab,Ttot_a,Ttot_b
        
an = (1.0/numpy.arange(1,sample_size)).sum()     
        
avg_TtotA = avg_TtotA*1.0/num_tree_replicates/an
avg_TtotB = avg_TtotB*1.0/num_tree_replicates/an

sys.stderr.write("Ttot: %g %g\n" % (avg_TtotA, avg_TtotB))