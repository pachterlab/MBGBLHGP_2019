#!/usr/bin/env python

import sys
from collections import Counter
import argparse
import pandas as pd

def makebug(t2g, trn, ecs, bus):
    def parseec(l):
        try: 
            tids = [int(i) for i in l.split(",")]
        except AttributeError: 
            tids = [l]
        genes = [trn.iloc[i]["gene"] for i in tids]
        genes = ",".join(map(str, list(Counter(genes).keys())))
        return genes

    print("Fixing up t2g..")
    t2g = t2g.set_index("trn")
    trn["gene"] = trn["trn"].map(t2g["gene"])
    print("Done.")

    print("Parsing matrix.ec..")
    ecs["gene"] = ecs.tid.map(parseec)
    print("Done.")

    print("Actually making the BUG file..")
    bus["gene"] = bus["ecs"].map(ecs["gene"])
    print("Done.")

    print("Bug file made internally")

    print("Cleaning up BUG file..")
    bug = bus.groupby(["bcs", "umi", "gene"])[["mul"]].sum()
    bug = bug.reset_index()
    print("Done.")
    return bug

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="converts a bus file to a bug file")
    parser.add_argument("--g", help="transcripts to genes file")
    parser.add_argument("--t", help="transcripts file")
    parser.add_argument("--e", help="equivalence classes file")
    parser.add_argument("--b", help="corrected sorted busfile")
    parser.add_argument("outdir", type=str, help="output directory for bug file")

    args = parser.parse_args()
    print("Loading files..")
    bus = pd.read_csv(args.b, header=None, names=["bcs" ,"umi", "ecs", "mul"], sep="\t")
    ecs = pd.read_csv(args.e, header=None, sep="\t", names=["ecs", "tid"])
    trn = pd.read_csv(args.t, header=None, names=["trn"])
    t2g = pd.read_csv(args.g, header=None, names=["trn", "gene", "gname"], sep="\t")
    print("files loaded. Making BUG file..")
    bug = makebug(t2g, trn, ecs, bus)
    print("BUG file made. Saving BUG file..")
    bug.to_csv(args.outdir + "bug.txt", sep="\t", header=False, index=False)
    print("BUG file saved!")
