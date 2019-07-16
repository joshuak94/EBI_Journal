#!/usr/bin/env python

import argparse
import random
import csv
from collections import defaultdict
import re

def sorted_nicely( l ):
    """ Sorts the given iterable in the way that is expected.

    Required arguments:
    l -- The iterable to be sorted.

    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Generate an input.tsv file for chipulate.")
    parser.add_argument("-b", "--bindingSites", type = int, default = 100, help = "The number of binding sites to be randomly distributed.")
    parser.add_argument("-c", "--chromosomes", type = str, nargs = '+', help = "The chromosomes to use in the distribution. Usage: -c 1 2 3 4 X", required = True)
    parser.add_argument("-d", "--directory", type = str, default = "./", help = "The output directory for input.tsv.")
    parser.add_argument("-l", "--lengthFile", type = str, help = "The full path to the chromosome lengths file.", required = True)
    parser.add_argument("-s", "--seed", action = "store_true", help = "Enable seed for reproducibility.")

    args = parser.parse_args()

    out_file = open("%s/input.tsv" % args.directory, "w+")
    out_file.write("chr\tstart\tend\tp_ext\tp_amp\tenergy_A\n")

    if args.seed:
        random.seed(1)

    chrLengths = {}
    results = defaultdict(list)
    with open(args.lengthFile) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            chrLengths[row[0]] = row[1]

    for i in range(args.bindingSites):
        # Get random chromosome
        chr = "chr%s" % random.choice(args.chromosomes)

        # Get random length of binding.
        bind = random.randint(150, 250)

        # Get random start on chromosome.
        start = random.randint(0, (int)(chrLengths[chr]) - bind - 200)

        # Get random extraction efficiency.
        p_ext = random.gauss(0.5, 0.1)

        # Get random p_amp
        p_amp = random.gauss(0.5, 0.1)

        # Get random energy_A
        energy_A = random.betavariate(2, 6)

        # Store in list.
        results[chr].append([start, bind, p_ext, p_amp, energy_A])

    for chr in sorted_nicely(results):
        prev = -1
        for i in sorted(results[chr]):
            # Write to output if unique.
            if i[0] != prev:
                out_file.write("%s\t%d\t%d\t%1.2f\t%1.2f\t%1.2f\n" % (chr, i[0], i[0] + i[1], i[2], i[3], i[4]))
            prev = i[0]
    out_file.close()
