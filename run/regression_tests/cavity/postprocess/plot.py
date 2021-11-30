#!/usr/bin/env python3
# coding: utf-8

import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from pathlib import Path

Re_number = 1000

# ghia
dfref = pd.read_csv(Path('../benchmark_data') / f'ghia_{Re_number}.csv')

# citylbm
result = Path('../build/run/result/cavity_V100_Nodes1_MPI4')

df_vx = []
for f in result.glob('**/output-station0.csv'):
    df_vx.append(pd.read_csv(f))
df_vx = pd.concat(df_vx)
df_vx = df_vx[df_vx['step'] == np.max(df_vx['step'])-1]
df_vx = df_vx.sort_values(['x'])

df_uy = []
for f in result.glob('**/output-station1.csv'):
    df_uy.append(pd.read_csv(f))
df_uy = pd.concat(df_uy)
df_uy = df_uy[df_uy['step'] == np.max(df_uy['step'])-1]
df_uy = df_uy.sort_values(['y'])

print(end='') # trick for jupyterlab


# plot
plt.rcParams['font.size'] = 16
#plt.rcParams['font.family'] = 'Times New Roman'
fig, ax = plt.subplots(figsize=(6,6))

# v
ax1 = ax.twinx()
ax1.plot(df_vx['x']+0.5, df_vx['v'], label='CityLBM', color='C0')
ax1.scatter(dfref['#ghia-x'], dfref['#ghia-v'], color='C1', label='Ghia, 1982')
ax1.set_xlim(0, 1)
ax1.set_ylim(-1, 1)
ax.set_xlabel('x')
ax1.set_ylabel('v')
ax.set_xticks(np.arange(0,1.5,0.5))
ax1.set_yticks(np.arange(-1,1.5,0.5))

# u
ax2 = ax.twiny()
ax2.plot(df_uy['u'], df_uy['y']+0.5, label='CityLBM', color='C0')
ax2.scatter(dfref['#ghia-u'], dfref['#ghia-y'], color='C1', label='Ghia, 1982')
ax2.set_xlabel('u')
ax.set_ylabel('y')
ax2.set_xlim(-1, 1)
ax2.set_ylim(0, 1)
ax.set_yticks(np.arange(0,1.5,0.5))
ax2.set_xticks(np.arange(-1,1.5,0.5))

ax2.legend()

plt.savefig(f'Re_{Re_number}.pdf', bbox_inches='tight')
