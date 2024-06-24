import os
import pandas as pd
import numpy as np
from brainsmash.mapgen.base import Base
from brainsmash.mapgen.eval import base_fit

dir = "/Users/isabellehales/Desktop"
os.chdir(dir)

# input and output path
path="/Users/isabellehales/Desktop/hoftmanlab/"

# function to compute n surrogates for a given map
# https://brainsmash.readthedocs.io/en/latest/gettingstarted.html#parcellated-surrogate-maps
def bransmash_generate(n:int, map:str, mat:str):
    print("generating "+str(n)+" permutations")
    print(map)
    print(mat)
    base = Base(x=map, D=mat)
    out = base(n=n)
    print("done")
    return out

AHBA_data = pd.read_excel("AHBA_500aparc_fs5.3_expmatrix_MSN_combat_cohend_HC_CHR_051623.xlsx")
msn_map = AHBA_data["cohend_HC_CHR"]
msn_labels = AHBA_data["ROI"]
#print(msn_labels)



m_map = msn_map.to_numpy()
#adding NaN to first index to account for unkownn region, excluded from AHBA data but in distance matrix
nm = np.insert(m_map, 0, "NaN")



# get permutations
perm = bransmash_generate(n=10000, map=nm, mat=path+"lh.500.aparc.distmat.txt")
# save as CSV
pd.DataFrame(perm).to_csv(path+"/PermTestAll.csv", header=False, index=False)




#m_map = m_map.reshape(-1,1)
#mf = m_map

#agene_d = pd.read_csv("gene_expression_n15043.csv", skiprows = [0], header = None, dtype = float)
#gene_d = agene_d.to_numpy()

#print(gene_d.shape)
