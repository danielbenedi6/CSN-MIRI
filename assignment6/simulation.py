import math
import random
from tqdm import tqdm
from multiprocessing import Process

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
    #for t in tqdm(range(1,tmax+1)):
    for t in range(1,tmax+1):
        new_stubs = random.sample(stubs, m0) + [n for _ in range(m0)]
        stubs += new_stubs
        n += 1
        
        for i in range(int(math.log10(tmax))):
            evolution[i][t] = evolution[i][t-1] + new_stubs.count(10**i+n0-1)

    return evolution,degree_sequence(stubs,n)

def random_attach_model(tmax, n0, m0):
    if n0 < 2:
        raise Exception("Argument n0 has to be larger than 2")
    
    deg_seq = [0 for _ in range(n0+tmax)]

    # Initialize with n0 nodes connected randomly
    # ensuring that all nodes have at least degree 1
    for n in range(n0):
        other_stub = random.randint(0,n0-2)
        if other_stub >= n:
            other_stub += 1
        deg_seq[n] += 1
        deg_seq[other_stub] += 1
    
    evolution = [[0 for _ in range(tmax+1)] for _ in range(int(math.log10(tmax)))]
    n = n0
    #for t in tqdm(range(1,tmax+1)):
    for t in range(1,tmax+1):
        deg_seq[n] += m0
        new_stubs = [ n for _ in range(m0)]
        for _ in range(m0):
            other_stub = random.randint(0,n-1)
            deg_seq[other_stub] += 1
            new_stubs += [other_stub]
        n += 1
        
        for i in range(int(math.log10(tmax))):
            evolution[i][t] = evolution[i][t-1] + new_stubs.count(10**i+n0-1)
    return evolution, deg_seq

def no_growth_model(tmax, n0, m0):
    stubs = []
    # Initialize with n0 nodes connected randomly
    # ensuring that all nodes have at least degree 1
    for n in range(n0):
        other_stub = random.randint(0,n0-2)
        if other_stub >= n:
            other_stub += 1
        stubs += [n, other_stub]

    evolution = [[0 for _ in range(tmax+1)] for _ in range(int(math.log10(n0)))]
    #for t in tqdm(range(1,tmax+1)):
    for t in range(1,tmax+1):
        n = random.randint(0,n0-1)
        new_stubs = random.sample(stubs, m0) + [n for _ in range(m0)]
        stubs += new_stubs
        
        for i in range(int(math.log10(n0))):
            evolution[i][t] = evolution[i][t-1] + new_stubs.count(10**i)
    
    return evolution,degree_sequence(stubs,n0)

def experiment(prefix, method, repetitions, tmax, n0, mo):
    for _ in range(rep):
        evolution,deg_seq = method(tmax, n0, m0)
        with open("./data/%s_deg_seq.txt"%(prefix),"a") as f:
            # Each line of the file represents one repetition of the experiment
            # and each line contains the degree sequence of the final graph
            for deg in deg_seq:
                f.write("%i "%(deg))
            f.write("\n")
        
        for i in range(len(evolution)):
            with open("./data/%s_evolution_%i.txt"%(prefix,10**i),"a") as f:
                # Each line of the file represents one repetition of the experiment
                # and each line contains the evolution of the i-th arrived node
                for deg in evolution[i]:
                    f.write("%i "%(deg))
                f.write("\n")

tmax = int(10e5)
n0 = 10
m0 = 3
rep = 10

if __name__ == '__main__':
    p1 = Process(target=experiment, args=("full_ba", ba_model, rep, tmax, n0, m0))
    p2 = Process(target=experiment, args=("random_attach", random_attach_model, rep, tmax, n0, m0))
    p3 = Process(target=experiment, args=("no_growth", no_growth_model, rep, tmax, n0**m0, m0))
    p1.start()
    p2.start()
    p3.start()
    p1.join()
    p2.join()
    p3.join()


