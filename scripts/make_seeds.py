import argparse
from itertools import combinations
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
    # parser.add_argument("-g", "--genomes", type=str, help="Paths to the genomes, comma delimited", required=True)
    parser.add_argument("-n", "--num-seeds", type=int, help="Number of spaced seeds to generate", required=True)
    parser.add_argument("-k", "--seed-length", type=int, help="Length of seeds", required=True)
    parser.add_argument("-w", "--weight", type=int, help="Number of weighted elements per seed", required=True)
    parser.add_argument("-e", "--entropy-bits", type=int, required=True,
                        help="Specify number of bits per seed for entropy calculation")
    parser.add_argument("-s", "--random-seed", type=int, help="Define a seed (default: 42)", default=42)
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


def calculate_entropy(seed, s_size):
    bit_counts = {}
    for i in itertools.product(['0', '1'], repeat=s_size):
        key = ''.join(i)
        bit_counts[key] = 0
    for i in range(len(seed) - s_size + 1):
        bit = seed[i:i+s_size]
        bit_counts[bit] += 1
    # print(bit_counts)
    # Calculate Shannon entropy
    entropy = 0
    for bit, count in bit_counts.items():
        if count == 0:
            continue
        num_substrings = len(seed) - s_size + 1
        entropy -= float(count)/float(num_substrings) * math.log(float(count)/float(num_substrings), 2)
    return entropy


def calculate_moment(seed, pivot):
    distances = []
    for idx, i in enumerate(seed):
        if i == '1':
            distances.append(abs(pivot - idx)**2)
    return sum(distances)


def nCr(n, r):
    return int(math.factorial(n)/math.factorial(r)/math.factorial(n-r))


def get_random_seeds(k, w, num, random_seed=False):
    if random_seed:
        np.random.seed(random_seed)
    if num > comb(k, w):
        raise ValueError('num cannot be greater than k choose w')
    # Choose the indices to be 1
    seeds = set()
    population = [i for i in range(k)]
    for _ in range(num):
        repeated_seed = True
        while repeated_seed:
            seed = ['']*k
            indices = np.random.choice(population, size=w, replace=False)
            for i in range(k):
                if i in indices:
                    seed[i] = '1'
                else:
                    seed[i] = '0'
            seed = "{}{}{}".format("1", "".join(seed), "1")
            repeated_seed = True if seed in seeds else False
        seeds.add(seed)
    return list(seeds)


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


def clean_header(i):
    return i[1:len(i)-1].replace('\'', '')


def main():
    args = parse_args()
    k = args.seed_length - 2
    w = args.weight - 2
    num_seeds = args.num_seeds
    seeds = get_random_seeds(k, w, num_seeds, random_seed=args.random_seed)
    calculate_entropy_vect = np.vectorize(calculate_entropy, excluded=['s_size'])
    entropies_1 = calculate_entropy_vect(seeds, 2)
    entropies_2 = calculate_entropy_vect(seeds, 3)
    entropies_3 = calculate_entropy_vect(seeds, 4)
    data = pd.DataFrame({'2bit': entropies_1, '3bit': entropies_2, '4bit': entropies_3},
                 index=seeds,
                 columns=['2bit', '3bit', '4bit'])
    data.to_csv(args.output, sep='\t')
    print(data.describe())

# def old_main():
#     args = parse_args()
#     genomes_paths = args.genomes.split(",")
#     genomes = {}
#     for genome_path in genomes_paths:
#         genomes[genome_path] = pyfaidx.Fasta(genome_path)[0][0:]
#     k = args.seed_length
#     w = args.weight
#     output = args.output
#     output_file = open(output, 'w+', buffering=1)
#     num_seeds = args.num_seeds
#
#     # seeds = get_random_seeds(k, w, 10)
#     print("Seeds generated")
#     calculate_entropy_vect = np.vectorize(calculate_entropy, excluded=['s_size'])
#     seeds = ['100'*20, '101000'*10, '100010'*10]
#     entropies = calculate_entropy_vect(seeds, args.entropy_bits)
#     # seeds = sample_seeds(seeds, entropies, num_seeds)
#     data_cols = list(combinations(genomes.keys(), 2))
#     entropies = calculate_entropy_vect(seeds, args.entropy_bits)
#     entropy_data = pd.DataFrame(entropies, index=seeds, columns=['entropies'])
#     print(entropy_data)
#     kmer_data = np.zeros(shape=(num_seeds, len(genomes)), dtype=object)
#     genome_kmers = {}
#     index = 0
#     print('seeds', '\t'.join([clean_header(str(i)) for i in combinations(genomes.keys(), 2)]), 'entropies', sep='\t', file=output_file)
#     for seed in seeds:
#         for genome_name, genome in genomes.items():
#             genome_kmers[genome_name] = determine_words(genome, seed)
#         pairwise_intersections = []
#         for comb in combinations(genome_kmers.values(), 2):
#             pairwise_intersections.append(str(len(set.intersection(*comb))))
#         print(seeds[index], '\t'.join(pairwise_intersections), str(entropies[index]), sep='\t', file=output_file)
#        index += 1

if __name__ == '__main__':
    main()
