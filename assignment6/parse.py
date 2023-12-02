import numpy as np
import os
import tqdm

def means(filename):
    f_in = open(filename, "r")
    f_out = open(filename+".mean","w")

    avg = np.fromstring(f_in.readline(), sep=" ")
    n = 1
    for line in f_in:
        avg += np.fromstring(line, sep=" ")
        n += 1
    avg = avg/n

    print(*avg, file=f_out)
    
def mean_spectrum(filename):
    f_in = open(filename, "r")
    f_out = open(filename+".mean","w")

    avg = np.zeros((1,))
    n = 0

    for line in f_in:
        k, f = np.unique(np.fromstring(line, sep=" ", dtype=int), return_counts=True)
        if max(k) >= len(avg):
            avg = np.pad(avg, max(k) - len(avg), "constant", constant_values=0)
        avg[k] += f
        n += 1
    avg = avg/n

    for i in range(len(avg)):
        if avg[i] > 0:
            print(f"{i} {avg[i]}", file=f_out)

for result in tqdm.tqdm(os.listdir("./data")):
    if ".mean" in result or ".csv" in result:
        continue
    #if "evolution" in result:
    #  means("./data/" + result)
    if "deg_seq" in result:
      print(result)
      mean_spectrum("./data/"+result)
