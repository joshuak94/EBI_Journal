#!/usr/bin/env python
import argparse
import subprocess
import sys
import os
from multiprocessing import Pool
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt

def runningMean(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]

def pearson(shift):
    p = subprocess.Popen(["wiggletools", "pearson", "strict", pos_file, "shiftPos", str(shift), neg_file],
                     stdout = subprocess.PIPE, stderr = subprocess.STDOUT)

    out, err = p.communicate()
    return shift, float(out.rstrip())

def pearson_masc(shift):
    command = "wiggletools pearson strict trimFill trim %s shiftPos %d %s unit %s trimFill trim %s shiftPos %d %s shiftPos %d unit %s" % (map_file, shift, map_file, pos_file, map_file, shift, map_file, shift, neg_file)
    pearson = float(subprocess.check_output(command.split()))
    return shift, pearson

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Calculate pearson correlation values for various shift lengths.")
    parser.add_argument("-c", "--cores", type = int, default = 1, help = "The number of cores to parallelize with.")
    parser.add_argument("-i", "--initial", type = int, default = 0, help = "The initial shift value to calculate.")
    parser.add_argument("-e", "--end", type = int, default = 400, help = "The final shift value to calculate.")
    parser.add_argument("-s", "--step", type = int, default = 1, help = "The step size for shift values.")
    parser.add_argument("-p", "--pos", type = str, help = "The name of the positive strand file.", required = True)
    parser.add_argument("-n", "--neg", type = str, help = "The name of the negative strand file.", required = True)
    parser.add_argument("-o", "--out", type = str, help = "The name of the file to save to.", required = True)
    parser.add_argument("-d", "--mapdir", type = str, help = "The relative path to the mappability file.", required = True)
    parser.add_argument("-r", "--readlen", type = int, help = "The read length for the bam files. Obtain via 'samtools view file.bam | head -1000 | awk '{print length($10)}' | sort -u'", required = True)
    parser.add_argument("-m", "--smooth", action = "store_true", help = "Enable smoothing.")
    args = parser.parse_args()
    global pos_file
    global neg_file
    global map_file
    pos_file = args.pos
    neg_file = args.neg
    map_file = args.mapdir
    out_file = open(args.out, "w+")
    out_file.write("Shift\tValue\n")
    pool = Pool(processes = args.cores)
    pearson_vals = np.array(pool.map(pearson_masc, range(args.initial, args.end + args.step, args.step)))
    maxPearson = np.nanmax(pearson_vals[:, 1])
    while (maxPearson == pearson_vals[0, 1] and pearson_vals.size > 2):
        pearson_vals = np.delete(pearson_vals, 0, 0)
        maxPearson = np.nanmax(pearson_vals[:, 1])
    if (pearson_vals.size == 2):
        sys.stderr.write("Error: Fragment length likely invalid. Continuously downward correlation trend.")
        out_file.close()
        sys.exit()
    if (args.smooth):
        print("Smoothing enabled")
        pearson_vals[:, 1] = runningMean(pearson_vals[:, 1], 15)
        maxPearson = np.nanmax(pearson_vals[:, 1])
    maxShift = pearson_vals[np.where(pearson_vals[:, 1] == maxPearson)[0][0], 0] + args.readlen
    print("Max pearson value: %2.6f\n" % maxPearson)
    print("Corresponding shift: %d\n" % maxShift)
    np.savetxt(out_file, np.vstack((pearson_vals[:, 0] + args.readlen, pearson_vals[:, 1])).T, fmt = "%d\t%2.6f", newline = "\n")
    if (args.smooth):
        plt.plot(pearson_vals[:-15, 0] + args.readlen, pearson_vals[:-15, 1], 'bo', markersize = 2, label = "Max shift: %d" % maxShift)
    else:
        plt.plot(pearson_vals[:, 0] + args.readlen, pearson_vals[:, 1], 'bo', markersize = 2, label = "Max shift: %d" % maxShift)
    plt.ylabel("Correlation values")
    plt.xlabel("Shift values")
    plt.legend()
    plt.savefig("%s.pdf" % args.out, format = "pdf")
    out_file.close()
