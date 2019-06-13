#!/usr/bin/env python

import argparse
import pandas as pd

def getumicelldist(bug):
   gb = bug.groupby("umi")["bcs"].nunique()
   return gb 

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="estimate N, number of umis from which we are sampling.")
    parser.add_argument("--b", help="bug file in text format")
    parser.add_argument("outdir", type=str, help="output directory for Nest.txt")
    args = parser.parse_args()
    
    print("Loading bug file..")
    bug = pd.read_csv(args.b, header=None, names=["bcs", "umi", "gene", "mul"], sep="\t")
    print("Done.")
    
    print("Getting distribution of number of cells per umi..")
    umicelldist = getumicelldist(bug)
    print("Saving..")
    umicelldist.to_csv(args.outdir + "umicelldist.txt", header=None, sep="\t")
