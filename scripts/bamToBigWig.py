#!/usr/bin/env python3
import argparse
import subprocess
import os
import shlex

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Submit wiggletools signal normalizer jobs for a series of files.")
    parser.add_argument("input", nargs = "+")

    args = parser.parse_args()

    for f in args.input:
        wig_file = f[:-4] + ".wig"
        bw_file = f[:-4] + ".bw"
        command = "bsub -M 8192 -R \"rusage[mem=8192]\" -o bamToBW.bjob.{}.out -e bamToBW.bjob.{}.err 'wiggletools write {}.wig {} && wigToBigWig {} http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes {} && rm {}'".format(f[:-4], f[:-4], wig_file, f, wig_file, bw_file, wig_file)
        subprocess.run(shlex.split(command, posix = False))
