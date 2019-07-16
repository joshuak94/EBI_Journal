#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os
import multiprocessing
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import itertools
import shlex
import math

################################################################################
# This is a reimplementation of align2rawsignal, originally written in matlab by Anshul Kundaje (https://github.com/akundaje/align2rawsignal).
# This reimplementation was written by Joshua Kim (joshua.kim@fu-berlin.de) in July 2019.
# It uses WiggleTools as the main tool to calculate various statistics given input files.
#
# USAGE:
# This tool should be called as follows: ./wiggletoolsSigNorm.py -i INPUT.bam -m MAPPABILITY.bw -c 1
# Here are some tips to running this tools as efficiently as possible:
# 1. The fragment length calculation can be run in parallel, but more cores (-c) must be specified
#    * If specifying more than 1 core, ensure enough memory is allocated. Generally, 1 Gb RAM per core.
# 2. The input *must* be in bam file format, and it is strongly encouraged that the mappability file be in bw file format.
# 3. The read length calculator looks at the sequence string of the INPUT file. If the sequence string is not present,
#    this will default to 1. Therefore, the option exists (-r) to manually give the read length and skip this calculation.
# 4. The assembly is assumed to be hg38/GRCh38 by default. If this is not the case, this can also be specified with -a.
#
# Note that this tool relies on a mappability wig/bw track for the specific read length of the given input file.
# This can be generated however you would like, but it is recommended to use GenMap (https://github.com/cpockrandt/genmap)
# by Chris Pockrandt, as that is the tool that was used during development.
################################################################################
def calcReadLen(input):
    command = "samtools view {}".format(input)
    p1 = subprocess.Popen(shlex.split(command), stdout = subprocess.PIPE)
    command = "head -1000"
    p2 = subprocess.Popen(shlex.split(command), stdin = p1.stdout, stdout = subprocess.PIPE)
    command = "awk '{print length($10)}'"
    p3 = subprocess.Popen(shlex.split(command), stdin = p2.stdout, stdout = subprocess.PIPE)
    command = "sort -ur"
    return int(subprocess.check_output(shlex.split(command), stdin = p3.stdout).split()[0])

def preprocess(args):
    url = "http://hgdownload.soe.ucsc.edu/goldenPath/{}/bigZips/{}.chrom.sizes".format(args.assembly, args.assembly)
    for strand in ("pos", "neg"):
        filename = "{}{}.bam".format(args.input[:-3], strand)
        splitter = "samtools view -h %s %s %s" % ("-F" if strand == "pos" else "-f", "20" if strand == "pos" else "16", args.input)
        sort = "samtools sort -o {} -@ {} -".format(filename, args.cores)
        convert = "{} -clip stdin {} {}.bw".format(args.wigToBigWig if args.wigToBigWig else "wigToBigWig", url, filename[:-4])
        # Delete file if it exists.
        if os.path.isfile(filename): os.remove(filename)
        if os.path.isfile(filename + ".bai"): os.remove(filename + ".bai")
        if os.path.isfile(filename[:-3] + "bw"): os.remove(filename[:-3] + "bw")
        # Split bam into pos and neg.
        try:
            if args.verbose: print("\tSplitting {} strand.".format(strand))
            p = subprocess.Popen(splitter.split(), stdout = subprocess.PIPE)
            if args.verbose: print("\tSorting bam file.")
            subprocess.run(shlex.split(sort), stdin = p.stdout)
            if args.verbose: print("\tGenerating index.")
            subprocess.run(["samtools", "index", filename], check = True)
            if args.verbose: print("\tGenerating bigwig file.")
            p = subprocess.Popen(["wiggletools", filename], stdout = subprocess.PIPE)
            subprocess.run(convert.split(), stdin = p.stdout, check = True)
        except subprocess.CalledProcessError as err:
            sys.exit(err)

def pearson_masc(args, shift):
    pos_file = "{}pos.bw".format(args.input[:-3])
    neg_file = "{}neg.bw".format(args.input[:-3])
    command = "wiggletools pearson strict trimFill trim %s shiftPos %d %s unit %s trimFill trim %s shiftPos %d %s shiftPos %d unit %s" % (args.map, shift, args.map, pos_file, args.map, shift, args.map, shift, neg_file)
    if args.naive: command = "wiggletools pearson strict %s shiftPos %d %s" % (pos_file, shift, neg_file)
    pearson = float(subprocess.check_output(shlex.split(command)))
    return pearson, shift

def runningMean(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]

def fragLengthCalc(args, readlen):
    pool = multiprocessing.Pool(processes = args.cores)
    pearson_vals = np.array(pool.starmap(pearson_masc, list(zip(itertools.repeat(args), range(args.begin, args.end + 1)))))
    pool.close()
    pool.join()
    maxPearson = np.nanmax(pearson_vals[:, 0])
    maxShift = pearson_vals[np.where(pearson_vals[:, 0] == maxPearson)[0][0], 1] + readlen
    if not args.naive:
        while (maxPearson == pearson_vals[0, 0] and pearson_vals.size - (2*args.smooth + 1) > 2):
            pearson_vals = np.delete(pearson_vals, 0, 0)
            maxPearson = np.nanmax(pearson_vals[:, 0])
        if (pearson_vals.size - (2*args.smooth + 1) == 2): sys.exit("Error: Fragment length likely invalid. Continuously downward correlation trend.\n")
        pearson_vals[:, 0] = runningMean(pearson_vals[:, 0], args.smooth)
        maxPearson = np.nanmax(pearson_vals[args.smooth:-args.smooth, 0])
        maxShift = pearson_vals[np.where(pearson_vals[args.smooth:-args.smooth, 0] == maxPearson)[0][0], 1] + readlen
    if args.verbose:
        out_file = open(args.out[:-3] if args.out else args.input[:-4] + "_fl.txt", "w+")
        out_file.write("Shift\tValue\n")
        if args.naive:
            plt.plot(pearson_vals[:, 1] + readlen, pearson_vals[:, 0], 'bo', markersize = 2, label = "Max shift: %d" % maxShift)
            np.savetxt(out_file, np.vstack((pearson_vals[:, 1] + readlen, pearson_vals[:, 0])).T, fmt = "%d\t%2.6f", newline = "\n")
        else:
            plt.plot(pearson_vals[args.smooth:-args.smooth, 1] + readlen, pearson_vals[args.smooth:-args.smooth, 0], 'bo', markersize = 2, label = "Max shift: %d" % maxShift)
            np.savetxt(out_file, np.vstack((pearson_vals[args.smooth:-args.smooth, 1] + readlen, pearson_vals[args.smooth:-args.smooth, 0])).T, fmt = "%d\t%2.6f", newline = "\n")
        plt.ylabel("Correlation values")
        plt.xlabel("Shift values")
        plt.legend()
        plt.savefig("%s.pdf" % args.out[:-3] if args.out else args.input[:-4], format = "pdf")
        out_file.close()
    return maxShift

def signalNorm(args, shift, readlen):
    pos_file = "{}pos.bw".format(args.input[:-3])
    neg_file = "{}neg.bw".format(args.input[:-3])
    wig_file = "{}.norm.wig".format(args.out[:-3] if args.out else args.input[:-4])
    url = "http://hgdownload.soe.ucsc.edu/goldenPath/{}/bigZips/{}.chrom.sizes".format(args.assembly, args.assembly)
    try:
        if args.verbose: print("\tCalculating positive and negative mapped read counts.")
        command = "samtools idxstats {}".format(pos_file[:-2] + "bam")
        p = subprocess.Popen(shlex.split(command), stdout = subprocess.PIPE)
        command = "awk 'BEGIN {sum = 0} {sum += $3} END {print sum}'"
        pos_reads = int(subprocess.check_output(shlex.split(command), stdin = p.stdout))
        command = "samtools idxstats {}".format(neg_file[:-2] + "bam")
        p = subprocess.Popen(shlex.split(command), stdout = subprocess.PIPE)
        command = "awk 'BEGIN {sum = 0} {sum += $3} END {print sum}'"
        neg_reads = int(subprocess.check_output(shlex.split(command), stdin = p.stdout))
        total_mapped_reads = pos_reads + neg_reads
        if args.verbose: print("\t\tPositive read count: {}.\n\t\tNegative read count: {}.".format(pos_reads, neg_reads))

        if args.verbose: print("\tCalculating total mappable positions.")
        command = "wiggletools AUC {}".format(args.map)
        mappable_genome_size = float(subprocess.check_output(shlex.split(command)))
        if args.verbose: print("\t\tMappable genome size: {}.".format(mappable_genome_size))
        total_ratio = (2*mappable_genome_size)/(total_mapped_reads*readlen)
        if args.verbose: print("\t\tTotal ratio: {}.".format(total_ratio))

        expected_file = "{}.expected.wig".format(args.out[:-3] if args.out else args.input[:-4])
        if args.verbose: print("\tCreating expected distribution file {}.".format(expected_file))
        command = "wiggletools write {} winsum {} sum shiftPos {} {} {}".format(expected_file, args.windowLength, shift - readlen, args.map, args.map)
        subprocess.run(shlex.split(command), check = True)

        if args.verbose: print("\tCalculating cutoff for regions of low mappability.")
        command = "wiggletools maxI {}".format(expected_file)
        threshold = (2 * args.windowLength) / 4
        minMaxLcmVal = 0.5 * 2 * args.windowLength

        if args.verbose: print("\tNormalizing signal into file {}".format(args.out if args.out else args.input[:-4] + ".norm.bw"))
        # command = "wiggletools write {} trim mult gt {} {} {} : ratio scale {} winsum {} sum {} shiftPos {} {} : max {} const {} {}".format(wig_file, threshold, expected_file, args.map, total_ratio, args.windowLength, pos_file, shift - readlen, neg_file, expected_file, minMaxLcmVal, expected_file)
        command = "wiggletools write {} ratio scale {} winsum {} sum {} shiftPos {} {} : max {} const {} {}".format(wig_file, total_ratio, args.windowLength, pos_file, shift - readlen, neg_file, expected_file, minMaxLcmVal, expected_file)
        p = subprocess.run(shlex.split(command), check = True)
        command = "{} -clip {} {} {}".format(args.wigToBigWig if args.wigToBigWig else "wigToBigWig", wig_file, url, args.out if args.out else args.input[:-4] + ".norm.bw")
        subprocess.run(shlex.split(command), check = True)
    except subprocess.CalledProcessError as err:
        print("Error: {}".format(err))

def cleanup(args):
    for strand in ("pos", "neg"):
        file_pre = "{}{}.".format(args.input[:-3], strand)
        for ext in ("bam", "bam.bai", "bw"):
            if args.verbose: print("\tRemoving {}".format(file_pre + ext))
            os.remove(file_pre + ext)
    expected_file = "{}.expected.wig".format(args.out[:-3] if args.out else args.input[:-4])
    wig_file = "{}.norm.wig".format(args.out[:-3] if args.out else args.input[:-4])
    os.remove(expected_file)
    os.remove(wig_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Create signal tracks for a ChIP-seq BAM file which have been normalized.\nRequires samtools and wigToBigWig to both exist in the PATH.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", type = str, required = True, help = "The input BAM file to use.")
    parser.add_argument("-m", "--map", type = str, help = "The relative path to the mappability file.", required = True)
    parser.add_argument("-c", "--cores", type = int, required = True, help = "The number of cores to parallelize with.")
    parser.add_argument("-o", "--out", type = str, help = "The name of the normalized output file. Defaults to [INPUT].norm.bw.")
    parser.add_argument("-b", "--begin", type = int, default = 0, help = "The initial shift value to calculate.")
    parser.add_argument("-e", "--end", type = int, default = 400, help = "The final shift value to calculate.")
    parser.add_argument("-s", "--smooth", type = int, default = 5, help = "The window to smooth the fragment length calculation by.")
    parser.add_argument("-w", "--wigToBigWig", type = str, help = "Path to wigToBigWig command. Only needs to be supplied if wigToBigWig is not in your path.")
    parser.add_argument("-a", "--assembly", type = str, help = "The assembly used (in UCSC notation). This is used to fetch the chromosome length file from http://hgdownload.soe.ucsc.edu/goldenPath/ASSEMBLY/bigZips/ASSEMBLY.chrom.sizes and defaults to hg38", default = "hg38")
    parser.add_argument("-f", "--fragLength", type = int, help = "Optional parameter if you want to skip fragment length calculation.")
    parser.add_argument("-r", "--readLength", type = int, help = "Optional parameter if you want to skip read length calculation.")
    parser.add_argument("-d", "--windowLength", type = int, default = 300, help = "Window length to use for signal aggregation.")
    parser.add_argument("-v", "--verbose", action = "store_true", help = "Print verbosely. Will also generate additional files for the fragment length calculation.")
    parser.add_argument("-n", "--nocleanup", action = "store_true", help = "Skip cleanup")
    parser.add_argument("-g", "--naive", action = "store_true", help = "Calculate only the fragment length using naive pearson correlation.")
    args = parser.parse_args()

    if args.input[-3:] != "bam": sys.exit("Input argument must be a single bam file.")
    if args.out and args.out[-2:] != "bw": sys.exit("If giving an output filename, it must end with the extension .bw")
    if args.map[-3:] == "wig": print("Warning: Fragment length calculations will be much slower with a wig file instead of a bigwig file.")
    readlen = args.readLength if args.readLength else calcReadLen(args.input)
    if args.verbose: print("Read length is {}.".format(readlen))
    if args.fragLength:
        frag_length = args.fragLength
    else:
        if args.verbose: print("Preprocessing input BAM file...")
        preprocess(args)
        frag_length = fragLengthCalc(args, readlen)
    if args.verbose: print("Fragment length is {}.".format(frag_length))
    if not args.naive:
        if args.verbose: print("Normalizing signal...")
        signalNorm(args, frag_length, readlen)
    if not args.nocleanup:
        if args.verbose: print("Cleaning up files...")
        cleanup(args)
