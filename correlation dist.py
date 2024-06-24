import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
dir = "/Users/isabellehales/Desktop"
os.chdir(dir)


# input and output path
path="/Users/isabellehales/Desktop/"

# load gene data
plsload = pd.read_excel("pls.scores.AHBA.geneexp.CHR_HC.MSN.combat.500aparc.5.3.cohend.model.simpls.comp1_122823.xlsx", header = None, skiprows = 1)
plsv = plsload.iloc[:,0]
#plsvn = plsv.values.flatten()


# load permutations
permd = pd.read_csv("PermTestWOUnknown.csv", header = None)

# run correlations
#correlations = np.corrcoef(plsload, permd)
total_corr = np.array([])

for i in range(10000):
    curr_spin = permd.iloc[i,:]
    #curr_spin_n = curr_spin.values.flatten()
    cur_cor = np.corrcoef(plsv, curr_spin)[0, 1]
    total_corr = np.append(total_corr, cur_cor)


#print(total_corr.shape)
csv_file_path = path+'allcor.csv'
# Use np.savetxt to save the array as CSV
np.savetxt(csv_file_path, total_corr, delimiter=',')

final_c = np.abs(total_corr)

# Plot the distribution of correlation values
plt.plot(final_c, color='black')
plt.title('Initial Distribution of Correlation Values')
plt.xlabel('Permutation Spin Number')
plt.ylabel('Correlation Coefficient')
plt.show()


plt.hist(total_corr, color='black', bins = 100)
plt.title('Initial Distribution of Correlation Values')
plt.xlabel('Permutation Spin Number')
plt.ylabel('Correlation Coefficient')
plt.show()

