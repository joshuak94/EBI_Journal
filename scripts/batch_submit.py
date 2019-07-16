#!/usr/bin/env python3
import argparse
import subprocess
import os
import shlex

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Submit wiggletools signal normalizer jobs for a series of files.")
    parser.add_argument("-w", "--wiggletools", required = True, help = "Path to wiggletoolsSigNorm.py file.")
    parser.add_argument("-m", "--mappability", required = True, help = "Path to mappability file with correct read lengths.")
    parser.add_argument("-o", "--output", required = True, help = "The output directory for all the final files.")
    parser.add_argument("input", nargs = "+")

    args = parser.parse_args()

    if not os.path.isfile(args.wiggletools): parser.error("Invalid path to wiggletoolsSigNorm.py.")
    if not os.path.isfile(args.mappability): parser.error("Invalid path to the mappability file.")
    if not os.path.isdir(args.output): parser.error("Invalid path to the output directory.")

    for f in args.input:
        command = "bsub -M 32768 -R \"rusage[mem=32768]\" -n 32 -o {}/bjob.{}.out -e {}/bjob.{}.err {} -c 32 -v -m {} -i {} -o {}/output.{}.out".format(args.output, f[:-4], args.output, f[:-4], args.wiggletools, args.mappability, f, args.output, f[:-4])
        subprocess.run(shlex.split(command, posix = False))
