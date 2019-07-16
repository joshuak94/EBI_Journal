#!/usr/bin/env python3
import argparse
import subprocess
import sys
import random
import multiprocessing
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib.colors as pltc
import numpy as np
import itertools
all_colors = [k for k,v in pltc.cnames.items()]

def pairwise_pearson(arguments):
    erb, bw_files = arguments
    bw1, bw2 = bw_files
    print("Calculating pearson for {} and {}".format(bw1, bw2))
    command = "wiggletools pearson trim {} {} trim {} {}".format(erb, bw1, erb, bw2)
    try:
        val = float(subprocess.check_output(command.split()))
    except subprocess.CalledProcessError as err:
        print("Error: {}".format(err))
        sys.exit(1)
    return bw1, bw2, val

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Given a list of bigWig files, compute the pairwise pearson correlation coefficients and print results to a tab delimited file.")
    parser.add_argument("-o", "--output", default = "output.pdf", help = "Name of the output file.")
    parser.add_argument("-c", "--cores", default = 1, type = int, help = "Number of cores to multiprocess.")
    parser.add_argument("-e", "--erb", required = True, type = str, help = "Path to the Ensemble Regulatory Build bigwig file.")
    parser.add_argument("input", nargs = "+")
    args = parser.parse_args()

    vals = np.identity(len(args.input))
    # Fill upper triangular part.
    pairs = [(args.input[i], args.input[j]) for i in range(len(args.input)) for j in range(i + 1, len(args.input))]
    pool = multiprocessing.Pool(processes = args.cores)
    results = np.array(pool.map(pairwise_pearson, list(zip(itertools.repeat(args.erb), pairs))))
    pool.close()
    pool.join()

    for item in results:
        vals[args.input.index(item[0])][args.input.index(item[1])] = item[2]
    print("Current pearson array:\n{}".format(vals))

    # Calculate eigen vectors and eigen values of above, using upper triangle.
    eigen_vals, eigen_vectors = np.linalg.eigh(vals, 'U')

    # Reduce pearson correlation matrix
    reduced = np.matmul(vals, eigen_vectors[len(eigen_vectors) - 2:].T)
    # print("Eigen values: {}\nEigen vectors: {}".format(eigen_vals, eigen_vectors))

    # Calculate PCA.
    # pca = PCA(n_components = 2)
    # vals_r = pca.fit(vals).transform(vals)
    # print(vals_r)

    # Plot PCA.
    plt.figure()
    colors = random.sample(all_colors, len(args.input))
    lw = 2

    plt.scatter(x = reduced[:, 0], y = reduced[:, 1], alpha = .8, lw = lw)
    # plt.scatter(x = vals_r[:, 0], y = vals_r[:, 1], alpha = .8, lw = lw)
    # for color, i, target_name in zip(colors, range(len(args.input)), args.input):
    #     plt.scatter(vals_r[i, 0], vals_r[i, 1], color=color, alpha=.8, lw=lw,
    #                 label=target_name)
    # plt.legend(loc='best', shadow=False, scatterpoints=1)
    plt.title('PCA of dataset')
    plt.savefig("{}.pdf".format(args.output), format = "pdf")
