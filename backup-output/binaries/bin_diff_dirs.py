#!/usr/bin/env python3

import array
import os
import sys

def floats_fromfile(fname, verbose):
    with open(fname, 'rb') as f:
        a = array.array('f')
        a.fromfile(f, int(os.path.getsize(fname)/4))
        if verbose:
            print(fname + ': ' + str(min(a)) + " -- " + str(max(a)))
        return a

def diff_file(fname_left, fname_right, verbose):
    left  = floats_fromfile(fname_left, verbose)
    right = floats_fromfile(fname_right, verbose)
    min_diff = 99999.9
    max_diff = 0.0
    if len(left) != len(right):
        print("invalid size of binaries: " + str(len(left)) + " -- " + str(len(right)))
        print("")
    else:
        for (l, r) in zip(left, right):
            min_diff = min(min_diff, abs(l-r))
            max_diff = max(max_diff, abs(l-r))
        if verbose:
            print("diff: " + str(min_diff) + " -- " + str(max_diff))
            print("")
    return (min_diff, max_diff)


import itertools 
import argparse

def main():
    # arg
    parser = argparse.ArgumentParser()
    parser.add_argument('dir_left' , help="name of left directory" )
    parser.add_argument('dir_right', help="name of right directory")
    parser.add_argument('-v', '--verbose', action="store_true", help="print verbose info")
    parser.add_argument('-s', '--startswith', default="",     help="string where the target filename starts with")
    parser.add_argument('-e', '--endswith'  , default=".dat", help="string where the target filename ends with")
    parser.add_argument('--exclude_lbm', action="store_true", help="not evalate lbm value f")
    parser.add_argument('--exclude_xyz', action="store_true", help="not evalate coordinate")
    args = parser.parse_args()

    # dir
    dir_left = args.dir_left
    dir_right = args.dir_right
    if not os.path.isdir(dir_left) and not os.path.isdir(dir_right):
        print("not such directory(ies): " + dir_left + " or " + dir_right)

    # do
    if args.verbose:
        print("dirs: " + dir_left + " -- " + dir_right)

    files_left =  [f.name for f in os.scandir(args.dir_left ) if f.name.startswith(args.startswith) and f.name.endswith(args.endswith)]
    files_right = [f.name for f in os.scandir(args.dir_right) if f.name.startswith(args.startswith) and f.name.endswith(args.endswith)]
    if files_left != files_right:
        print("directories are not same")
        print(" continueing anyway")

    # files: files_left v files_right
    files = []
    files.extend(files_left )
    files.extend(files_right)
    files = list(set(files))
    sort(files)

    # exclude
    if args.exclude_lbm:
        files = [f for f in files if not f.startswith('lbm_')]
    if args.exclude_xyz:
        files = [f for f in files if not f.startswith('x_')]
        files = [f for f in files if not f.startswith('y_')]
        files = [f for f in files if not f.startswith('z_')]

    print("num of files = " + str(len(files)))

    # calc minmax foreach file
    min_diff = 99999.9
    max_diff = 0.0
    for fname in files:
        print(fname)
        ffname_left  = dir_left  + "/" + fname
        ffname_right = dir_right + "/" + fname
        tmp_min_diff,tmp_max_diff = diff_file(ffname_left, ffname_right, args.verbose)
        min_diff = min(min_diff, tmp_min_diff)
        max_diff = max(max_diff, tmp_max_diff)
        #print(fname_left + ": " + str(min_diff) + " -- " + str(max_diff))

    print("merged diff range: " + str(min_diff) + " -- " + str(max_diff))
    
    # final
    if max_diff < 1.e-8:
        return 0
    else:
        return 1
  
if __name__ == '__main__':
    sys.exit(main())

