#!/usr/bin/env python3
import argparse
import subprocess
import os
import shlex
import multiprocessing

def calculateStats(file):
    command = "wiggletools AUC CVI varI minI maxI {}".format(file)

    print("Calculating for {}.".format(file))
    auc, cvi, var, min, max = subprocess.check_output(shlex.split(command)).split()

    return file, auc, cvi, var, min, max

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Submit wiggletools signal normalizer jobs for a series of files.")
    parser.add_argument("-c", "--cores", type = int, required = True, help = "Number of cores to use")
    parser.add_argument("-o", "--out", type = str, required = True, help = "Output file name")
    parser.add_argument("input", nargs = "+")

    args = parser.parse_args()

    out_file = open(args.out, "w+")
    out_file.write("File name\tAUC\tCOV\tVariance\tMin\tMax\n")

    pool = multiprocessing.Pool(processes = args.cores)
    final = pool.map(calculateStats, args.input)

    for file, auc, cvi, var, min, max in final:
        out_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(file, float(auc), float(cvi), float(var), float(min), float(max)))
