#!usr/bin/env python3

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Normalize the signal of a track based on mappability.")
    parser.add_argument("-m", "--mappability", type = str, required = True, help = "The path to the mappability file, preferrably wig or bw.")
    parser.add_argument("-p", "--posbam", type = str, required = True, help = "The bam file containing only reads mapping to the positive strand.")
    parser.add_argument("-n", "--negbam", type = str, required = True, help = "The bam file containing only reads mapping to the negative strand.")
    parser.add_argument("-")
