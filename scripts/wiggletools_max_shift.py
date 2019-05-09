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

pos_file = ""
neg_file = ""

def pearson(shift):
    print("Processing %d" % shift)
    p = subprocess.Popen(["wiggletools", "pearson", pos_file, "shiftPos", str(shift), neg_file],
                         stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    out, err = p.communicate()
    return shift, float(out.rstrip())

# For wigCorrelate
def wigCorrelate(shift):
    print("Processing %d" % shift)
    tmpfile = "/scratch/joshjayk/shift%d.wig" % shift
    command = "wiggletools write %s shiftPos %s %s;../../wigCorrelate %s %s;rm %s" % (tmpfile, str(shift), neg_file, pos_file, tmpfile, tmpfile)
    p = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell = True)
    out, err = p.communicate()
    return shift, float(out.split("\t")[2].rstrip())

# The new method
def accurate(shift):
    print("Processing %d" % shift)
    tmpfile = "/scratch/joshjayk/shift%d.wig" % shift
    map_shift = 36 - shift
    command = "wiggletools write %s overlaps ../mappability/filtered_mus_musculus_all_chrom.map.wig shiftPos %s ../mappability/filtered_mus_musculus_all_chrom.map.wig;wiggletools pearson overlaps %s shiftPos %s %s overlaps %s %s" % (tmpfile, str(map_shift), tmpfile, str(shift), neg_file, tmpfile, pos_file)
    p = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, shell = True)
    out, err = p.communicate()
    return shift, float(out.rstrip())

def pearson_masc(shift):
    print("Processing %d" % shift)
    p = subprocess.Popen(["wiggletools", "pearson", pos_file, "shiftPos", str(shift), neg_file],
                         stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    out, err = p.communicate()
    return shift, float(out.rstrip())
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Calculate pearson correlation values for various shift lengths.")
    parser.add_argument("-c", "--cores", type = int, default = 1, help = "The number of cores to parallelize with.")
    parser.add_argument("-i", "--initial", type = int, default = 1, help = "The initial shift value to calculate.")
    parser.add_argument("-e", "--end", type = int, default = 400, help = "The final shift value to calculate.")
    parser.add_argument("-s", "--step", type = int, default = 1, help = "The step size for shift values.")
    parser.add_argument("-m", "--masc", action = "store_true", help = "Use MaSC formula for calculating neg strand coverage.")
    parser.add_argument("-w", "--wigCorrelate", action = "store_true", help = "Use wigCorrelate to calculate pearson values.")
    parser.add_argument("-a", "--accurate", action = "store_true", help = "Temporary to test the new method.")
    parser.add_argument("-p", "--pos", type = str, help = "The name of the positive strand file.", required = True)
    parser.add_argument("-n", "--neg", type = str, help = "The name of the negative strand file.", required = True)
    parser.add_argument("-o", "--out", type = str, help = "The name of the file to save to.", required = True)
    args = parser.parse_args()
    pos_file = args.pos
    neg_file = args.neg
    out_file = open(args.out, "w+")
    out_file.write("Shift\tValue\n")
    pool = Pool(processes = args.cores)
    if args.masc:
        print("MaSC set.")
        pearson_vals = np.array(pool.map(pearson_masc, range(args.initial, args.end + args.step, args.step)))
    elif args.wigCorrelate:
        print("wigCorrelate set.")
        pearson_vals = np.array(pool.map(wigCorrelate, range(args.initial, args.end + args.step, args.step)))
    elif args.accurate:
        print("New method.")
        pearson_vals = np.array(pool.map(accurate, range(args.initial, args.end + args.step, args.step)))
    else:
        print("Wiggletools set.")
        pearson_vals = np.array(pool.map(pearson, range(args.initial, args.end + args.step, args.step)))
    maxPearson = np.amax(pearson_vals, axis = 0)[1]
    maxShift = pearson_vals[np.where(pearson_vals[:, 1] == maxPearson), 0]
    print("Max pearson value: %2.6f\n" % maxPearson)
    print("Corresponding shift: %d\n" % maxShift)
    np.savetxt(out_file, pearson_vals, fmt = "%d\t%2.6f", newline = "\n")
    plt.plot(pearson_vals[:, 0], pearson_vals[:, 1], 'bo')
    plt.ylabel("Correlation values")
    plt.xlabel("Shift values")
    plt.annotate("Max shift: %d" % maxShift, xy = (maxShift, maxPearson), xytext = (maxShift, maxPearson + 0.1),
                  arrowprops = dict(facecolor = "black", shrink = 0.1))
    plt.savefig("%s.pdf" % args.out, format = "pdf")
    out_file.close()
