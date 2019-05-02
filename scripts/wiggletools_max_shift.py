import argparse
import subprocess
import sys
from multiprocessing import Pool

pos_file = ""
neg_file = ""

def pearson(shift):
    print("Processing %d" % shift)
    p = subprocess.Popen(["wiggletools", "pearson", pos_file, "shiftPos", str(shift), neg_file],
                         stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    out, err = p.communicate()
    return float(out.rstrip())

def pearson_masc(shift):
    print("Processing %d" % shift)
    p = subprocess.Popen(["wiggletools", "pearson", pos_file, "shiftPos", str(shift), neg_file],
                         stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    out, err = p.communicate()
    return float(out.rstrip())
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Calculate pearson correlation values for various shift lengths.")
    parser.add_argument("-c", "--cores", type = int, default = 1, help = "The number of cores to parallelize with.")
    parser.add_argument("-i", "--initial", type = int, default = 1, help = "The initial shift value to calculate.")
    parser.add_argument("-e", "--end", type = int, default = 400, help = "The final shift value to calculate.")
    parser.add_argument("-s", "--step", type = int, default = 1, help = "The step size for shift values.")
    parser.add_argument("-m", "--masc", action = "store_true", help = "Use MaSC formula for calculating neg strand coverage.")
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
        pearson_vals = pool.map(pearson_masc, range(args.initial, args.end + args.step, args.step))
    else:
        print("MaSC not set.")
        pearson_vals = pool.map(pearson, range(args.initial, args.end + args.step, args.step))
    print("Max pearson value: %2.6f\n" % (max(pearson_vals)))
    print("Corresponding shift: %d\n" % (pearson_vals.index(max(pearson_vals)) + 1))
    for i, val in enumerate(pearson_vals, start = 1):
        out_file.write("%d\t%2.6f\n" % (i, val))
    out_file.close()
