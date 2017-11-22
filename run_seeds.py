import argparse
from matplotlib_venn import venn3
from matplotlib_venn._venn3 import compute_venn3_subsets
from matplotlib import pyplot as plt
import designSS
import pyfaidx
import itertools
import math
import numpy as np
import pandas as pd
from scipy.special import comb


# def parse_args():
#     parser = argparse.ArgumentParser(
#         description="""Runs a number of spaced seeds against specified genomes
#                     to determine if seed entropy affects specificity""")
#     parser.add_argument("-i", "--input", type=str, help="Barcode whitelist path", required=True)
#     parser.add_argument("-m", "--molecules", type=str, help="Molecules/amplicons fasta file", required=True)
#     parser.add_argument("-o", "--output", type=str, help="Output file (default: stdout)")
#     parser.add_argument("-s", "--random-seed", type=int, help="Define a seed (default: no seed)")
#     args = parser.parse_args()
#     return args


def search_sequence(sequence, spaced_seed, kmer):
    weights_to_kmers = {}
    j = 0
    for i in range(len(spaced_seed)):
        if spaced_seed[i] == '1':
            weights_to_kmers[i] = kmer[j]
            j += 1

    if not len(weights_to_kmers) == len(kmer):
        raise ValueError('The length of the kmer must be the same as the weight of the spaced seed')
    for i in range(len(sequence) - len(spaced_seed)):
        seed_in_sequence = True
        for index, nc in weights_to_kmers.items():
            if not sequence[i + index] == nc:
                seed_in_sequence = False
                break
        if seed_in_sequence:
            return True
    return False


def determine_words(genome_sequence, spaced_seed):
    kmers = set()
    k = len(spaced_seed)
    weighted_indexes = []
    for i in range(k):
        if spaced_seed[i] == '1':
            weighted_indexes.append(i)
    for i in range(len(genome_sequence) - k + 1):
        kmer = ''
        for index in weighted_indexes:
            # Note: genome_sequence is a pyfaidx.Sequence
            kmer += str(genome_sequence[i:i+k][index])
        kmers.add(kmer)
    return kmers


def calculate_entropy(seed, s_size):
    smer_counts = {}
    for i in itertools.product(['0', '1'], repeat=s_size):
        key = ''.join(i)
        smer_counts[key] = 0
    for i in range(len(seed) - s_size + 1):
        smer = seed[i:i+s_size]
        smer_counts[smer] += 1
    # print(smer_counts)
    # Calculate Shannon entropy
    entropy = 0
    for smer, count in smer_counts.items():
        if count == 0:
            continue
        num_substrings = len(seed) - s_size + 1
        entropy -= float(count)/float(num_substrings) * math.log(float(count)/float(num_substrings), 2)
    return entropy


def get_random_seeds(k, w, num, random_seed=False):
    if random_seed:
        np.random.seed(random_seed)
    if num > comb(k, w):
        raise ValueError('num cannot be greater than k choose w')
    # Choose the indices to be 1
    seeds = []
    for _ in range(num):
        repeated_seed = True
        while repeated_seed:
            seed = ''
            indices = np.random.choice([i for i in range(k)], size=w, replace=False)
            for i in range(k):
                if i in indices:
                    seed += '1'
                else:
                    seed += '0'
            repeated_seed = True if seed in seeds else False
        seeds.append(seed)
    return seeds


def get_random_kmers(k, num, random_seed=False):
    bases = ['A', 'C', 'G', 'T']
    if random_seed:
        np.random.seed(random_seed)
    kmers = []
    for _ in range(num):
        repeated_kmer = True
        while repeated_kmer:
            kmer = "".join(np.random.choice(bases, size=k, replace=True))
            repeated_kmer = True if kmer in kmers else False
        kmers.append(kmer)
    return kmers


def main():
    # seed = designSS.design_seed()
    bases = ['A', 'C', 'G', 'T']
    ecoli_k12 = pyfaidx.Fasta('e_coli.fa.gz')[0][0:]
    ecoli_BW25113 = pyfaidx.Fasta('e_ecoli_BW25113.fa')[0][0:]
    yeast = pyfaidx.Fasta('Pichia_sorbitophila.fa')[0][0:]
    genomes = {'E. Coli k12': ecoli_k12, "E. Coli BW25113": ecoli_BW25113, "Pichia sorbitophila": yeast}
    # print(search_sequence('AAACAAAAGTAACG', '1000011', 'CGT'))
    k = 60
    w = 20
    num_seeds = 100

    # kmers = get_random_kmers(w, 20)
    seeds = get_random_seeds(k, w, num_seeds)
    calculate_entropy_vect = np.vectorize(calculate_entropy, excluded=['s_size'])
    entropies = calculate_entropy_vect(seeds, 3)

    data = pd.DataFrame(entropies, index=seeds, columns=['entropies'])
    genome_kmers = {}
    for seed in seeds:
        for genome_name, genome in genomes.items():
            genome_kmers[genome_name] = determine_words(genome, seed)
        venn3(genome_kmers.values(), genome_kmers.keys())
        plt.savefig("seed_{}.png".format(seed))
    # kmer_values = list(genome_kmers.values())
    # compute_venn3_subsets(kmer_values[0], kmer_values[1], kmer_values[2])
    # `data.to_csv('genome_{}.tsv'.format(genome_name), sep='\t')


if __name__ == '__main__':
    main()
