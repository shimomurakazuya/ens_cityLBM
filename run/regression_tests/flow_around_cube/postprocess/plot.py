#!/usr/bin/env python3
# coding: utf-8

# Plot results of FlowAroundCube

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter
import pandas as pd
import numpy as np
import scipy as sp
import seaborn as sns
from pathlib import Path

#### set path:
result = Path('../build/run/result/flow_cube_V100_Nodes10_MPI40/io/output')
benchmark = Path('../benchmark_data')

# wind xz plane
def wxz(u='u', dt=0.01, ds=1, steps=[100,160], icsv=1):
    fig, ax = plt.subplots(dpi=100)
    
    # cal
    df = []
    for f in result.glob(f'**/output-station{icsv}.csv'):
        df.append(pd.read_csv(f))
    df = pd.concat(df)
    df = df.query(f'step >= {steps[0]} and step <= {steps[1]}')
    dfg = df.groupby(['z', 'y', 'x'], as_index=False)
    dfa = dfg.mean()
    dfs = dfg.std()
    
    for x in [-.15, -.05, 0, .065, .1, .15, .2]:
        dfas = dfa[abs(dfa['x'] - x) < 0.02]
        z = dfas['z']
        ua = dfas[u]
        dfss = dfs[abs(dfs['x'] - x) < 0.02]
        us =  dfss[u]*ds
        ax_cal, = ax.plot(x + ua*dt, z, 'r-', marker=None, zorder=102)
        ax.fill_betweenx(z, x + (ua-us)*dt, x + (ua+us)*dt, color='r', alpha=0.2)
        ax.plot([x, x], [0, 1], 'k:')
    
    # ref
    dfref_xy = pd.read_csv(benchmark / 'velocity/TableW2-2.csv')
    ax_ref = ax.scatter(dfref_xy[u.upper()+' (m/s)']*dt + dfref_xy['X (mm)'] / 1000, dfref_xy['Z (mm)'] / 1000, )
    try:
        ax_ref = ax.errorbar(dfref_xy[u.upper()+' (m/s)']*dt*ds + dfref_xy['X (mm)'] / 1000, dfref_xy['Z (mm)'] / 1000, xerr=dfref_xy['sigma ' + u + '(m/s)']*dt, fmt='oC0', ecolor='gray')
    except KeyError:
        ax_ref = ax.errorbar(dfref_xy[u.upper()+' (m/s)']*dt*ds + dfref_xy['X (mm)'] / 1000, dfref_xy['Z (mm)'] / 1000, xerr=dfref_xy['sigma ' + u + ' (m/s)']*dt, fmt='oC0', ecolor='gray')
    
    # finalize
    ax.add_patch(plt.Rectangle(xy=[-0.05,0], width=0.1, height=0.1, facecolor='gray', edgecolor='k', zorder=99))
    ax.set_xlabel(fr'X [m], X + {dt}{u.upper()} [m]')
    ax.set_ylabel('Z [m]')
    ax.set_ylim(0, 0.5)
    ax.set_xlim(-0.2, 0.25)
    ax.set_title("Vertical wind profile on spanwise cross-section (Y=0), " + u.upper())
    ax.legend([ax_ref, ax_cal], ['Exp', 'Calc'], bbox_to_anchor=(1.01, 1))
    plt.savefig(f'wind_xz_{u}.pdf', bbox_inches='tight')
    plt.close(fig)

wxz('u', dt=0.02, ds=1)
wxz('v', dt=0.02, ds=1)
wxz('w', dt=0.02, ds=1)

# wind xy plane
def wxy(u='u', dt=0.01, icsv=0, steps=[100, 160]):
    fig, ax = plt.subplots(dpi=100)
    
    # cal
    df = []
    for f in result.glob(f'**/output-station{icsv}.csv'):
        df.append(pd.read_csv(f))
    df = pd.concat(df)
    df = df.query(f'step >= {steps[0]} and step <= {steps[1]}')
    dfg = df.groupby(['z', 'y', 'x'], as_index=False)
    dfa = dfg.mean()
    dfs = dfg.std()
    
    for x in [-.05, 0, .065, .1, .15, .2, .25]:
        dfas = dfa[abs(dfa['x'] - x + 0.005) < 0.01]
        y = dfas['y']
        ua = dfas[u]
        dfss = dfs[abs(dfs['x'] - x + 0.005) < 0.01]
        us = dfss[u]
        ax_cal, = ax.plot(x + ua*dt, y, 'r-', zorder=102)
        ax.fill_betweenx(y, x + dt*(ua-us), x + dt*(ua+us), color='r', alpha=0.2)
        ax.plot([x,x], [-1,1], 'k:')
    
    # ref
    dfref_xy = pd.read_csv(benchmark / 'velocity/TableW2-1.csv')
    ax_ref = ax.scatter(dfref_xy[u.upper() + ' (m/s)']*dt + dfref_xy['X (mm)'] / 1000, dfref_xy['Y (mm)'] / 1000)
    try:
        ax_ref = ax.errorbar(dfref_xy[u.upper()+' (m/s)']*dt + dfref_xy['X (mm)'] / 1000, dfref_xy['Y (mm)'] / 1000, xerr=dfref_xy['sigma ' + u + '(m/s)']*dt, fmt='oC0', ecolor='k')
    except KeyError:
        ax_ref = ax.errorbar(dfref_xy[u.upper()+' (m/s)']*dt + dfref_xy['X (mm)'] / 1000, dfref_xy['Y (mm)'] / 1000, xerr=dfref_xy['sigma ' + u + ' (m/s)']*dt, fmt='oC0', ecolor='k')

    # finalize
    ax.add_patch(plt.Rectangle(xy=[-0.05, -0.05], width=0.1, height=0.1, facecolor='gray', edgecolor='black', zorder=99))
    ax.set_xlabel(fr'X [m], X + {dt}{u.upper()} [m]')
    ax.set_ylabel('Y [m]')
    ax.set_xlim(-0.1, 0.3)
    ax.set_ylim(-0.2, 0.2)
    ax.set_title('Spanwise wind profile on horizontal cross-section (Z = H/2), ' + u.upper())
    ax.legend([ax_ref, ax_cal], ['Exp', 'Calc'], bbox_to_anchor=(1.01, 1))
    plt.savefig(f'wind_xy_{u}.pdf', bbox_inches='tight')
    plt.close(fig)
    
wxy('u', dt=0.02)
wxy('v', dt=0.02)
wxy('w', dt=0.02)


# scalar xy plane
def sxy(scalar='scalar0', xmax=2, icsv=3, steps=[100, 160]):
    fig, ax = plt.subplots(dpi=100)
    
    # cal
    df = []
    for f in result.glob(f'**/output-station{icsv}.csv'):
        df.append(pd.read_csv(f))
    df = pd.concat(df)
    df = df.query(f'step >= {steps[0]} and step <= {steps[1]} and x <= {xmax}')
    df = df.drop_duplicates(['y', 'x'])
    dfg = df.groupby(['y', 'x'], as_index=False)
    dfa = dfg.mean()
    dfs = dfg.std()
    
    sns.set_style('white')
    levels = [pow(10,i*0.2) for i in range(21)]
    ax.contour(df.pivot('y', 'x', scalar), levels=levels, colors=['black'], alpha=0.5, zorder=102)
    contour = ax.contour(df.pivot('y', 'x', scalar), levels=[1, 10, 100], colors=['red'], zorder=102)
    contour.clabel(fmt='%.0f', fontsize=12)
    contourf = ax.contourf(df.pivot('y', 'x', scalar), levels=levels, cmap='Blues', norm=LogNorm(), zorder=101)
    fig.colorbar(contourf, ticks=levels, ax=ax)
    
    # finalize
    ax.set_title('Contour of scalar concentration on the ground level, calc')
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    xsize = df.drop_duplicates(subset=['x']).shape[0]
    ysize = df.drop_duplicates(subset=['y']).shape[0]
    ax.set_xticks([0, xsize/2, xsize])
    ax.set_xticklabels(['0', '1', '2'], rotation='horizontal')
    ax.set_yticks([0, ysize/2, ysize])
    ax.set_yticklabels(['-0.5', '0', '0.5'])
    ax.add_patch(plt.Rectangle(xy=[-xsize/40, ysize/2-ysize/20], width=xsize/20, height=ysize/10, facecolor='gray', edgecolor='black', zorder=999))
    plt.savefig(f'scalar_contour_cal.pdf', bbox_inches='tight')
    plt.close(fig)
    
sxy()

def sxy_ref():
    fig, ax = plt.subplots(dpi=100)
    
    # ref

    def contour_df(df):
        df_x = df.columns[1:]
        xx = pd.to_numeric(df_x) / 1000
        df_y = df['Y(mm)'] / 1000
        df_val = df[df_x]
        levels = [pow(10,i*0.2) for i in range(21)]
        ax.contour(xx, df_y, df_val, levels=levels, colors=['black'], alpha=0.5, zorder=102)
        contour = ax.contour(xx, df_y, df_val, levels=[1, 10, 100], colors=['red'], zorder=102)
        contour.clabel(fmt='%.0f', fontsize=12)
        contourf = ax.contourf(xx, df_y, df_val, levels=levels, cmap='Blues', norm=LogNorm(), zorder=101)
        fig.colorbar(contourf, ticks=levels, ax=ax)
    df = pd.read_csv(benchmark / 'scalar/TableC3-1.csv')
    df2 = df[df['Y(mm)'] % 40 == 0]
    contour_df(df2)
    
    # finalize
    ax.set_title('Contour of scalar concentration on the ground level, exp')
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_xlim(0,2)
    ax.set_ylim(-0.5, 0.5)
    ax.set_xticks([0,1,2])
    ax.set_yticks([-0.5, 0, 0.5])
    ax.add_patch(plt.Rectangle(xy=[-0.05, -0.05], width=0.1, height=0.1, facecolor='gray', edgecolor='black'))
    
    plt.savefig(f'scalar_contour_exp.pdf', bbox_inches='tight')
    plt.close(fig)
    
sxy_ref()

def stat(scalar='scalar0', icsv=5, steps=[100,160]):
    # cal
    df = []
    for f in result.glob(f'**/output-station{icsv}.csv'):
        df.append(pd.read_csv(f))
    df = pd.concat(df)
    df = df.query(f'step >= {steps[0]} and step <= {steps[1]}')
    df = df.drop_duplicates(['y', 'x'])
    dfg = df.groupby(['y', 'x'], as_index=False)
    dfc = dfg.mean()

    dft = pd.read_csv(benchmark / 'scalar/TableC3-1.csv')
    dft_x = dft.columns[1:]
    xx = pd.to_numeric(dft_x) / 1000
    dft_y = dft['Y(mm)']
    rxx = []
    exx = []
    cxx = []
    cnt_cnan = 0
    for y in dft_y:
        dy = dft[dft['Y(mm)'] == y]
        y_real = y/1000
        for x in dft.columns[1:]:
            x_real = float(x)/1000

            dx = dy[x]
            d = None
            for dd in dx:
                d = dd
            
            calc = dfc
            calc = calc[abs(calc['y'] - y_real) < 0.001]
            calc = calc[abs(calc['x'] - x_real) < 0.001]
            calc = calc[scalar]
            
            c = None
            for cc in calc:
                c = cc
            
            if d is not None and not np.isnan(d):
                if c is not None and not np.isnan(c):
                    exx.append(d)
                    cxx.append(c)
                    r = d/c
                    rxx.append(r)
                else:
                    cnt_cnan = cnt_cnan + 1
    
    return np.array(cxx), np.array(exx), np.array(rxx), cnt_cnan
    
def qq():
    fig, ax = plt.subplots(dpi=100)
    
    cxx, exx, _, __ = stat()
    ax.scatter(exx, cxx, marker='.', label='exp. vs. calc.')
    tt = np.array([1e-2, 400])
    ax.plot(tt, tt, color='k', ls=':')
    ax.set_aspect('equal')
    ax.fill_between(tt, tt/2, tt*2, color='k', alpha=0.1, label='FAC2')
    ax.set_xlabel('experiment [ppm]')
    ax.set_ylabel('calculation [ppm]')
    ax.set_title('Scatter plot of scalar concentrations')
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_xlim(*tt)
    ax.set_ylim(*tt)
    ax.legend()
    plt.savefig(f'scalar_scatter.pdf', bbox_inches='tight')
    plt.close(fig)
    
qq()

def factor2():
    _, __, rxx, cnt_cnan = stat()
    a = rxx.size
    b = np.count_nonzero(( (rxx >= 0.5) & (rxx <= 2) ))
    print('factor2 score: {} ({}/{}) at cnt_calc_nan={}'.format(b/a, b, a, cnt_cnan))
    with open('fac2.txt', 'w') as f:
        print('factor2 score: {} ({}/{}) at cnt_calc_nan={}'.format(b/a, b, a, cnt_cnan), file=f)
    
factor2()
