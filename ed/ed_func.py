import numpy as np

def new_basis_num(sz,i,j):
    """
    Given sz array of spin states sz and indices i and j of spin sites to flip,
    return the new integer corresponding to the spin state after the flip
    """
    sztemp = np.copy(sz) #make shallow copy so sz array isnt modified
    if sztemp[i] < 0:
        sztemp[i] = 1/2
    else:
        sztemp[i] = -1/2
    if sztemp[j] < 0:
        sztemp[j] = 1/2
    else:
        sztemp[j] = -1/2
    new_num = 0 
    sztemp = np.flip(sztemp,0) #note the flip to restore natural ordering!
    for i in range(len(sztemp)):
        new_num += (2**i) * (abs(sztemp[i]-1/2)) #note the abs and -1/2
    return int(new_num)

def get_spin_val(basis_num):
    """
    given the current basis_num, what is the spin value at the current site
    """
    if basis_num%2 == 0:
        return 1/2
    else:
        return -1/2
        

def hamil_basis_num(x,dimer_sites,sz_0_sym=False):
    """
    x is the basis num, dimer_sites is a list of tuples containing coordinates of dimer sites
    sz_0_sym is if you intend to utilize sz=0 symmetry
    """
    
    sz = np.zeros(nn)
    pbc_y = nn - nx #number of sites difference for PBC in y
    basis_num = x
    
    for i in range(nn): #loop over number of sites
        sz[i] = get_spin_val(basis_num) # map to spin up and down
        basis_num = basis_num//2 # '//' to denote integer division
    sz = np.flip(sz,0) #flip to go with natural basis ordering

    if sz_0_sym:
        if sum(sz) != 0: 
            return None 
    
    data = []
    row = []
    col = []
    
    basis_num = x
    
    for i in range(nn): 
        
        #do the x-direction bonds first
        if (i+1)%nx != 0: #if it is not the last site of every row
            
            # -------------------- diagonal elements --------------------
            data.append(J1*sz[i]*sz[i+1]) 
            row.append(basis_num)
            col.append(basis_num)
        
            # ----------------- off diagonal elements --------------------
            if sz[i] != sz[i+1]: #sz on neighboring sites must be different for hopping element to be non-zero
                basis_num_aft_hop = new_basis_num(sz,i,i+1)
                data.append(J1/2) #action of J1/2(S+S- + S-S+) on S=1/2 gives eigenvalue of J1/2
                row.append(basis_num_aft_hop)
                col.append(basis_num)
        
        else: #implement pbc in x direction
            if periodic_x:
                # -------------------- diagonal elements --------------------
                data.append(J1*sz[i]*sz[i-nx+1])  
                row.append(basis_num)
                col.append(basis_num)

                # ----------------- off diagonal elements --------------------
                if sz[i] != sz[i-nx+1]: 
                    basis_num_aft_hop = new_basis_num(sz,i,i-nx+1)
                    data.append(J1/2) 
                    row.append(basis_num_aft_hop)
                    col.append(basis_num)
            
        #do the y-direction bonds now
        if i < nn-nx: #if its not the last site of every column
            
            # -------------------- diagonal elements --------------------
            data.append(J1*sz[i]*sz[i+nx])
            row.append(basis_num)
            col.append(basis_num)
            
            # ----------------- off diagonal elements --------------------
            if sz[i] != sz[i+nx]:
                basis_num_aft_hop = new_basis_num(sz,i,i+nx)
                data.append(J1/2) 
                row.append(basis_num_aft_hop)
                col.append(basis_num)
                
        else: #if it is last site of every column, then implement pbc
            if periodic_y:
                # -------------------- diagonal elements --------------------
                data.append(J1*sz[i]*sz[i-pbc_y])
                row.append(basis_num)
                col.append(basis_num)
                
                #off diagonal elements
                if sz[i] != sz[i-pbc_y]:
                    basis_num_aft_hop = new_basis_num(sz,i,i-pbc_y)
                    data.append(J1/2) 
                    row.append(basis_num_aft_hop)
                    col.append(basis_num)
                
    #interdimer coupling
    for (i,j) in dimer_sites:
        
        #diagonal
        data.append(J2 *sz[i] * sz[j])
        row.append(basis_num)
        col.append(basis_num)
        
        #off diagonal
        if sz[i] != sz[j]:
            basis_num_aft_hop = new_basis_num(sz,i,j)
            data.append(J2/2)
            row.append(basis_num_aft_hop)
            col.append(basis_num)
    
    return np.array(data), np.array(row), np.array(col),np.array(basis_num)

