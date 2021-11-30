#! /usr/bin/env python3

import numpy as np
import pandas as pd

def make_csv():
    # unit: meter
    def make_csv_sub(i, xr, yr, zr, xrc=[], yrc=[], zrc=[], *, ddx=0):
        xr = np.array(xr)
        yr = np.array(yr)
        zr = np.array(zr)
        xrc = np.array(xrc)
        yrc = np.array(yrc)
        zrc = np.array(zrc)
        size = xr.size * yr.size * zr.size + xrc.size * yrc.size * zrc.size
        print('input size:', size)
        xyz = np.zeros((size, 3))
        cnt = 0
        for z in zr:
            for y in yr:
                for x in xr:
                    xyz[cnt, 0] = x + ddx
                    xyz[cnt, 1] = y + ddx
                    xyz[cnt, 2] = z + ddx
                    cnt = cnt + 1
        for z in zrc:
            for y in yrc:
                for x in xrc:
                    xyz[cnt, 0] = x + ddx
                    xyz[cnt, 1] = y + ddx
                    xyz[cnt, 2] = z + ddx
                    cnt = cnt + 1
        np.savetxt('input-{i}.csv'.format(i=i), xyz, fmt='%f', delimiter=',')
    
    # 
    #
    dx = 5e-3 # 5mm
    dxc = dx * 4
    xxx = np.array([-0.15,  -0.05, 0    , 0.065, 0.1  , 0.15 , 0.2  , 0.25 ])
    xxx = xxx - dx/2 * ((xxx + 1e-30) / (abs(xxx) + 1e-30))
    print(xxx)
    # 0. wind xy plane
    make_csv_sub(0, 
        # fine
        xxx,
        np.arange(-0.2 - dx/2, 0.2 + dx/2 + dx, dx),
        [0.05],
        # coarse
        xxx,
        np.arange(-0.7 - dxc/2, 0.7 + dxc/2 + dxc, dxc),
        [0.05]
    )

    # 1. wind xz plane
    make_csv_sub(1,
        # fine
        xxx,
        [0],
        np.arange(0 + dx/2, 0.5 + dx/2, dx),
        # coarse
        xxx,
        [0],
        np.arange(0 + dxc/2, 1.72 + dxc/2, dxc)
        )

    # inflow yz planes
    make_csv_sub(2,
        # fine
        [-10, -5, -0.15],
        np.arange(-0.5 - dx/2, 0.5 + dx*3/2, 10*dx),
        np.arange(0 - dx/2, 0.5 + dx*3/2, 10*dx),
        # coarse
        [-10, -5, -0.15],
        np.arange(-0.72 - dxc/2, 0.72 + dxc*3/2),
        np.arange(0 - dxc/2, 1.72 + dxc/2 + dxc, 5*dxc)
        )

    # scalar xy plane (ground level)
    make_csv_sub(3,
        np.arange(0 + dx/2, 3 + dx*3/2, 10*dx),
        np.arange(-0.5 - dx/2, 0.5 + dx*3/2, 10*dx),
        [dx/2],
        )

    # scalar xy plane (ground level), for Q-Q plot
    make_csv_sub(5,
        1e-3*np.array([65, 100, 150, 200, 250, 300, 450, 680, 980, 1280, 1580, 2180, 2780]), # ond_data/scalar/TableC3-1.csv
        1e-3*np.arange(-480, 480, 20),
        [dx/2],
        )

    # scalar xz plane
    make_csv_sub(4,
        np.arange(-dx/2, 3 + dx*3/2, 10*dx*3/2),
        [0],
        np.arange(dx/2, 0.5 + dx*3/2, 10*dx*3/2),
        )

if __name__ == '__main__':
    make_csv()
