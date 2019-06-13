#!/usr/bin/env python

import sys
import argparse
import pandas as pd
import numpy as np
from scipy.optimize import fsolve

def solve(n, *data): # dont know n or k, know x # of umis that collide bw genes, suppI_g support of I for the gene, suppI_c support of I for cell
    x, suppI_g, suppI_c = data
    s = 0
    
    for i in suppI_g:
        s += i/(n-i)
    
    return suppI_c - (n-suppI_c)*s - x

def estNk(ubug):
    Nest = []
    umipercell = ubug.groupby("bcs")[["umi"]].count()
    umipercellpergene = ubug.groupby(["bcs", "gene"])[["umi"]].count()
    counter = 0
    for cellwewant in umipercell.umi.nlargest(50).keys():
        if counter%10==0:
            print(counter)
        
        cellbug = ubug[ubug["bcs"] == cellwewant]
        suppI_g  = cellbug.groupby(["bcs", "gene"])["umi"].nunique().values
        suppI_c = cellbug["umi"].nunique()
        
        x = cellbug.groupby(["bcs", "umi"])["gene"].nunique()
        x = x[x>1].shape[0]

        sol = fsolve(solve, 250000, args=(x, suppI_g, suppI_c))
        
        Nest.append(int(sol[0]))
        counter += 1
    return Nest

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="estimate N, number of umis from which we are sampling.")
    parser.add_argument("--b", help="bug file in text format")
    parser.add_argument("outdir", type=str, help="output directory for Nest.txt")
    args = parser.parse_args()

    print("Loading bug file..")
    bug = pd.read_csv(args.b, header=None, names=["bcs", "umi", "gene", "mul"], sep="\t")
    ubug = bug[bug["gene"].map(lambda l: "," not in str(l))]
    
    Nest = estNk(ubug)

    Nest = ",".join(map(str, Nest))

    with open(args.outdir + "Nest.txt", "w") as f:
        f.write("{}".format(Nest))
