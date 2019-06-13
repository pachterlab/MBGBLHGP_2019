#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np


def estkg(N, suppI_g):
    return np.log(1-suppI_g/N)/np.log(1-1/N)

def getlossratio(bug, N):
    ubug = bug[bug["gene"].map(lambda l: "," not in str(l))]
    suppI_g = pd.DataFrame(ubug.groupby(["bcs", "gene"])["umi"].nunique())
    suppI_g["kg"] = suppI_g["umi"].map(lambda x: estkg(N, x))
    suppI_g["loss"] = suppI_g["kg"]-suppI_g["umi"]
    suppI_g["loss_ratio"] = suppI_g["loss"]/suppI_g["kg"]
    return suppI_g["loss_ratio"]

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Compute the loss ratio per gene for all genes in a cell across all cells")
    parser.add_argument("--b", help="bug file in text format")
    parser.add_argument("--N", help="Nest.txt")
    parser.add_argument("outdir", type=str, help="output directory for umispercell.txt")
    args = parser.parse_args()
    
    print("Loading bug file..")
    bug = pd.read_csv(args.b, header=None, names=["bcs", "umi", "gene", "mul"], sep="\t")
    print("Done.")
    with open(args.N, 'r') as f:
        Nest = np.array(list(map(int, f.readlines()[0].split(","))))
    
    print("Counting UMIs....")
    lossratios = getlossratio(bug, Nest.mean())
    print("Saving..")
    with open(args.outdir + "lossratio.txt", "w") as f:
        f.write("{}".format(",".join(map(str, list(lossratios.values)))))
