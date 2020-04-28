"""
example driver code to run a QMC simulation of the dimer model
"""
import numpy as np
import sys

#sys.path.append('/home/users/ntu/aheng004/dimer/python/src')
sys.path.append('/path/to/src')

from mc_sse_dimer import *

if __name__ == "__main__":

    param = np.loadtxt('params.txt')
    if len(param) != 5:
        sys.exit("insufficient parameters")

    j1 = float(param[0])
    j2 = float(param[1])
    h = float(param[2])
    beta = float(param[3])
    L = int(param[4])


    test = mc_sse_dimer(j1,j2,h,beta,L)
    test.main()
    