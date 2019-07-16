#!/usr/bin/env python3

import argparse
import subprocess
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Preprocess a sorted and indexed BAM file.")
    parser.add_argument("-i", "--input", type = str, help = "The input bam file to process.", required = True)
    parser.add_argument("-m", "--mappability", type = str, help = "Mappability file to use.", required = True)
    parser.add_argument("-b", "--wigToBigWig", type = str, help = "Full path name to wigToBigWig executable", required = True)
    parser.add_argument("-o", "--overlaps", action = "store_true", help = "Use overlaps instead of trimFill.")
    parser.add_argument("-c", "--lengths", type = str, help = "Path to the chrom length file.", required = True)
    args = parser.parse_args()

    for strand in ("pos", "neg"):
        bam = open("%s%s.bam" % (args.input[:-3], strand), "wb")
        wt_command = "overlaps" if args.overlaps else "trimFill"
        wig_file = "%s%s.wig" % (bam.name[:-3], wt_command)
        bw_file = "%sbw" % (wig_file[:-3])
        wt = "wiggletools write %s %s %s %s" % (wig_file, wt_command, args.mappability, bam.name)
        splitter = "samtools view -bh %s %s %s" % ("-F" if strand == "pos" else "-f", "20" if strand == "pos" else "16", args.input)
        # Split bam into pos and neg.
        try:
            print("Processing %s." % bam.name)
            print("\tCommand run: %s" % splitter)
            subprocess.run(splitter.split(), check = True, stdout = bam)
            print("Generating index.")
            print("\tCommand run: samtools index %s" % bam.name)
            subprocess.run(["samtools", "index", bam.name], check = True)
        except subprocess.CalledProcessError as err:
            sys.exit(err)

        # Run wiggletools overlaps or trimFill.
        try:
            print("Running wiggletools %s for %s." % (wt_command, bam.name))
            print("\tCommand run: %s" % wt)
            subprocess.run(wt.split(), check = True)
        except subprocess.CalledProcessError as err:
            sys.exit(err)

        # Convert to bigwig.
        try:
            print("Converting %s to bw file." % wig_file)
            print("\tCommand run: %s %s %s %s" % (args.wigToBigWig, wig_file, args.lengths, bw_file))
            subprocess.run([args.wigToBigWig, wig_file, args.lengths, bw_file], check = True)
        except subprocess.CalledProcessError as err:
            sys.exit(err)

        bam.close()
