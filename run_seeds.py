import argparse
from itertools import combinations
from matplotlib_venn import venn3
from matplotlib import pyplot as plt
import designSS
import pyfaidx
import itertools
import math
import numpy as np
import pandas as pd
from scipy.special import comb


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Runs a number of spaced seeds against specified genomes
                    to determine if seed entropy affects specificity.""")
    parser.add_argument("-g", "--genomes", type=str, help="Paths to the genomes, comma delimited", required=True)
    parser.add_argument("-n", "--num-seeds", type=int, help="Number of spaced seeds to generate", required=True)
    parser.add_argument("-k", "--seed-length", type=int, help="Length of seeds", required=True)
    parser.add_argument("-w", "--weight", type=int, help="Number of weighted elements per seed", required=True)
    parser.add_argument("-e", "--entropy-bits", type=int, required=True,
                        help="Specify number of bits per seed for entropy calculation")
    parser.add_argument("-s", "--random-seed", type=int, help="Define a seed (default: no seed)")
    parser.add_argument("-o", "--output", type=str, help="Name of output .tsv file", required=True)
    args = parser.parse_args()
    return args


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
    w = len(weighted_indexes)
    for i in range(len(genome_sequence) - k + 1):
        word = np.zeros(shape=w, dtype=str)
        for j, weighted_index in zip(range(w), weighted_indexes):
            # Note: genome_sequence is a pyfaidx.Sequence
            word[j] = str(genome_sequence[i:i+k][weighted_index])
        kmers.add("".join(word))
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


def sample_seeds(seeds, entropies, sample_size):
    select_indexes = np.logspace(0, math.log10(len(entropies)-1), sample_size).astype(int)
    sorted_indexes = np.argsort(entropies)
    indexes = np.sort(sorted_indexes[select_indexes])
    return np.array(seeds)[indexes]


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
    args = parse_args()
    genomes_paths = args.genomes.split(",")
    genomes = {}
    for genome_path in genomes_paths:
        genomes[genome_path] = pyfaidx.Fasta(genome_path)[0][0:]
    k = args.seed_length
    w = args.weight
    output = args.output
    num_seeds = args.num_seeds

    seeds = get_random_seeds(k, w, 1000)
    print("Seeds generated")
    calculate_entropy_vect = np.vectorize(calculate_entropy, excluded=['s_size'])
    entropies = calculate_entropy_vect(seeds, args.entropy_bits)
    seeds = sample_seeds(seeds, entropies, num_seeds)
    data_cols = list(combinations(genomes.keys(), 2))
    entropies = calculate_entropy_vect(seeds, args.entropy_bits)
    entropy_data = pd.DataFrame(entropies, index=seeds, columns=['entropies'])
    print(entropy_data)
    kmer_data = np.zeros(shape=(num_seeds, len(genomes)), dtype=object)
    genome_kmers = {}
    index = 0
    for seed in seeds:
        for genome_name, genome in genomes.items():
            genome_kmers[genome_name] = determine_words(genome, seed)
        pairwise_intersections = []
        for comb in combinations(genome_kmers.values(), 2):
            pairwise_intersections.append(len(set.intersection(*comb)))
        print(pairwise_intersections)
        kmer_data[index] = np.array(pairwise_intersections)
        index += 1
        kmer_data_df = pd.DataFrame(kmer_data, index=seeds, columns=data_cols)
        combined_data = kmer_data_df.join(entropy_data)
        combined_data.to_csv(output, sep='\t')


if __name__ == '__main__':
    main()
