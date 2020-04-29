import time
from scipy import linalg
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
from joblib import Parallel, delayed
import sys
from pathlib import Path
import os

import ss_ed_func.py

if __name__ == "__main__":

    t_start=time.time()

    # Parameters to define for ED #######################
    periodic_x = True
    periodic_y = True
    dimer_sites = [(0,6),(2,8),(7,11),(9,13),(10,16),(12,18),(17,1),(19,3)] #manually encode the dimer sites
    unit_cells = 2
    nx = 2*unit_cells if periodic_x else 2*unit_cells+1
    ny = 2*unit_cells if periodic_y else 2*unit_cells+1
    nn = nx * ny
    J1 = 1
    J2 = 4 #ensures the dimer singlet ground state
    basis_max = 2**nn
    #####################################################################

    path = os.getcwd()
    output_folder = "/output"
    if not os.path.exists(path+output_folder):
        os.makedirs(path+output_folder)
    os.chdir(path+output_folder)

    out = open("output.txt",'w')
    result_file = open("energy.txt",'w')

    result = Parallel(n_jobs=-1,verbose=1,pre_dispatch='2*n_jobs')(delayed(hamil_basis_num)(i,dimer_sites) for i in range(basis_max))

    matrix_ele = []; rows = []; cols = []; basis_list = []
    for i in range(len(result)):
        matrix_ele.extend(result[i][0])
        rows.extend(result[i][1])
        cols.extend(result[i][2])
        basis_list.append(result[i][3])
    del result

    unique = np.unique(rows)
    new_indices = np.arange(0,len(unique),1)
    arr = np.zeros(max(rows) + 1, dtype=new_indices.dtype)
    arr[unique] = new_indices
    rows = arr[rows]
    arr = np.zeros(max(cols) + 1, dtype=new_indices.dtype)
    arr[unique] = new_indices
    cols = arr[cols]
    size = len(np.unique(rows))
    H = csr_matrix((matrix_ele,(rows,cols)),shape=(size,size),dtype=complex)
    t_eig1 = time.time()
    eigval, eigvec = eigsh(H,k=20,which='SA')
    t_eig2 = time.time()
    out.write("Time for diagonalizing unique matrix is {0:.3f}s \n".format(t_eig2-t_eig1))
    out.flush()

    np.savez_compressed("eigvals",a=eigval)
    np.savez_compressed("eigvecs",a=eigvec)
    t_end = time.time()
    out.write("Total time of running script {0:.3f}s".format(t_end-t_start))
    out.close()