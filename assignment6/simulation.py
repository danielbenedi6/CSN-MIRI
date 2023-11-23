import math
import random
from tqdm import tqdm

def degree_sequence(stubs, n):
    deg_seq = [0 for _ in range(n)]
    for stub in stubs:
        deg_seq[stub] += 1
    return deg_seq

def ba_model(tmax, n0, m0):
    if n0 < 2:
        raise Exception("Argument n0 has to be larger than 2")
    
    stubs = []
    # Initialize with n0 nodes connected randomly
    # ensuring that all nodes have at least degree 1
    for n in range(n0):
        other_stub = random.randint(0,n0-2)
        if other_stub >= n:
            other_stub += 1
        stubs += [n, other_stub]
    
    evolution = [[0 for _ in range(tmax+1)] for _ in range(int(math.log10(tmax)))]
    n = n0
    for t in tqdm(range(1,tmax+1)):
        new_stubs = random.sample(stubs, m0) + [n for _ in range(m0)]
        stubs += new_stubs
        n += 1
        
        for i in range(int(math.log10(tmax))):
            evolution[i][t] = evolution[i][t-1] + new_stubs.count(10**i+n-1)

    return evolution,degree_sequence(stubs,n)

tmax = int(10e5)
rep = 100

for _ in range(rep):
	evolution,deg_seq = ba_model(tmax, 10, 3)
	with open("./data/full_ba_deg_seq.txt","a") as f:
	    for deg in deg_seq:
	        f.write("%i "%(deg))
	    f.write("\n")
	
	for i in range(int(math.log10(tmax))):
	    with open("./data/full_ba_evolution_%i.txt"%(10**i),"a") as f:
	        for deg in evolution[i]:
	            f.write("%i "%(deg))
	        f.write("\n")

