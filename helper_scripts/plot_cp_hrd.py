#!/usr/bin/env python3

import pandas as pd
import re
import sys
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

## Usage:
## python /mnt/data/Novaseq/.software/plot_cp_hrd.py {input.alt_sol} {input.lpp} {input.hrd} {output.pdf}

## Reads alternative entries table from sequenza
alt = pd.read_csv(sys.argv[1],sep=" ")
entries = alt.sort_values(by='SLPP',ascending=False)

## Reads lpp values created by the cp_heatmap.R script
lpp = pd.read_csv(sys.argv[2],sep=" ")


## Reads hrd values created by the cp_heatmap.R script
dat = []
dat = pd.read_csv(sys.argv[3],sep=" ",header=None,usecols=[1,2,3,4,5,6])
data = pd.DataFrame(dat)
data.columns = ["LOH","TAI","LST","HRD","Cellularity","Ploidy"] 
data = data[data["Cellularity"] > 0]

mat = pd.pivot_table(data,values='HRD',columns='Ploidy',index='Cellularity')
mat = mat.iloc[::-1]


## Init figure 
fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(111)
plt.title(sys.argv[1])


## Actual heatmap 
heat = sns.heatmap(mat,cmap='coolwarm',vmin=10,vmax=80,center=42,annot=True,square=False,cbar=False,ax=ax, fmt='g')



## Set stars to actual sequenza solutions
if len(entries) > 0:
    i = 0
    for _,ent in entries.iterrows():

        if ent["SLPP"] != ent["SLPP"]: # if SLPP is NaN
            continue

        if i == 0:
            ax.scatter(((ent["ploidy"])/0.1)-10+.5, ((1-ent["cellularity"])/0.05)+0.07, marker='*', s=1200, color='yellow')
            ax.text(((ent["ploidy"])/0.1)-10+.15, ((1-ent["cellularity"])/0.05)+0.07, str(int(ent["HRDsum"])),size=12)
        else:
            ax.scatter(((ent["ploidy"])/0.1)-10+.5, ((1-ent["cellularity"])/0.05)+0.07, marker='*', s=900, color='darkgrey')
            ax.text(((ent["ploidy"])/0.1)-10+.15, ((1-ent["cellularity"])/0.05)+0.07, str(int(ent["HRDsum"])),size=12)
        i += 1

## Print the contours for the lpp values
np.seterr(divide='ignore')
lpp = np.log10(lpp)
lpp = -1*lpp
lpp = lpp.iloc[::-1]
lpp = lpp.T

X = list(lpp.columns)
X = [10*(x-1) for x in X]
Y = list(lpp.index)
Y = [21*(1-float(y)) for y in Y]
Z = lpp.values.tolist()

plt.contour(X,Y,Z,40,cmap='gist_gray',alpha=0.7)


## export final image

heat.get_figure().savefig(sys.argv[4],dpi=50,bbox_inches='tight')
