#!/usr/bin/env python3

import argparse
import subprocess
import sys
import shlex

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Preprocess a sorted and indexed paired end ChIP seq file.")
    parser.add_argument("-i", "--input", type = str, help = "The input bam file to process.", required = True)

    args = parser.parse_args()

    unpaired_bam = open("unpaired_%sbam" % args.input[:-3], "wb")
    open_command = "samtools view -h %s" % args.input
    unpair_command = "awk '{OFS=\"\t\"}{gsub(\"147|1171|1107|81|145|83\", \"16\", $2); gsub(\"163|1187|1123|65|129|99\", \"0\", $2); print $0}'"
    write_command = "samtools view -hSb -"

    # Unpair the bam file.
    print("Open bam file")
    p = subprocess.Popen(open_command.split(), stdout = subprocess.PIPE)
    print("Pipe to awk")
    p2 = subprocess.Popen(shlex.split(unpair_command), stdin = p.stdout, stdout = subprocess.PIPE)
    print("Write to unpaired.bam")
    subprocess.run(write_command.split(), stdin = p2.stdout, stdout = unpaired_bam)
    unpaired_bam.close()

    # # Split the bam file.
    # print("Splitting bam file")
    # pos_command = "samtools view -bh -F 20 %s" % unpaired_bam.name
    # neg_command = "samtools view -bh -f 16 %s" % unpaired_bam.name
    # pos_bam = open("%spos.bam" % args.input[:-3], "wb")
    # neg_bam = open("%sneg.bam" % args.input[:-3], "wb")
    #
    # p = subprocess.Popen(pos_command.split(), stdout = pos_bam)
    # p2 = subprocess.Popen(neg_command.split(), stdout = neg_bam)
    #
    # p.wait()
    # p2.wait()
    #
    # pos_bam.close()
    # neg_bam.close()
    #
    # # Sort and ndex.
    # print("Sorting")
    # p = subprocess.Popen(["samtools", "sort", "-o", "pos.sorted.%s" % args.input, pos_bam.name])
    # p = subprocess.Popen(["samtools", "sort", "-o", "neg.sorted.%s" % args.input, neg_bam.name])
    # p.wait()
    # p2.wait()
    # print("Indexing pos and neg bam files")
    # subprocess.Popen(["samtools", "index", "pos.sorted.%s" % args.input])
    # subprocess.Popen(["samtools", "index", "neg.sorted.%s" % args.input])
