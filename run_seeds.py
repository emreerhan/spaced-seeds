import designSS
import pyfaidx
import itertools
import math
import numpy as np
import pandas as pd
from scipy.special import comb


def search_sequence(sequence, spaced_seed, kmer):
    weights_to_kmers = {}
    j = 0
    for i in range(len(spaced_seed)):
        if spaced_seed[i] == "1":
            weights_to_kmers[i] = kmer[j]
            j += 1

    if not len(weights_to_kmers) == len(kmer):
        raise ValueError("The length of the kmer must be the same as the weight of the spaced seed")
    for i in range(len(sequence) - len(spaced_seed)):
        seed_in_sequence = True
        for index, nc in weights_to_kmers.items():
            if not sequence[i + index] == nc:
                seed_in_sequence = False
                break
        if seed_in_sequence:
            return True
    return False


def determine_spaced_kmers(genome_sequence, spaced_seeds):
    kmers_array = []
    for spaced_seed in spaced_seeds:
        kmers = set()
        k = len(spaced_seed)
        weighted_indexes = []
        for i in range(k):
            if spaced_seed[i] == "1":
                weighted_indexes.append(i)
        for i in range(len(genome_sequence) - k):
            kmer = ""
            for index in weighted_indexes:
                # Note: genome_sequence is a pyfaidx.Sequence
                kmer += str(genome_sequence[i:i+k][index])
            kmers.add(kmer)
        kmers_array.append(kmers)
    return kmers_array


def calculate_entropy(seed, s_size):
    smer_counts = {}
    for i in itertools.product(["0", "1"], repeat=s_size):
        key = "".join(i)
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
        raise ValueError("num cannot be greater than k choose w")
    # Choose the indices to be 1
    seeds = []
    for _ in range(num):
        repeated_seed = True
        while repeated_seed:
            seed = ""
            indices = np.random.choice([i for i in range(k)], size=w, replace=False)
            for i in range(k):
                if i in indices:
                    seed += "1"
                else:
                    seed += "0"
            repeated_seed = True if seed in seeds else False
        seeds.append(seed)
    return seeds


def main():
    # seed = designSS.design_seed()
    bases = ['A', 'C', 'G', 'T']
    ecoli_k12 = pyfaidx.Fasta('e_coli.fa.gz')[0][0:]
    ecoli_BW25113 = pyfaidx.Fasta('e_ecoli_BW25113.fa')[0][0:]
    yeast = pyfaidx.Fasta('Pichia_sorbitophila.fa')[0][0:]
    # print(search_sequence("AAACAAAAGTAACG", "1000011", "CGT"))
    k = 6
    w = 2
    num_seeds = 1
    kmers = ["AC", "CT", "GA"]
    seeds = get_random_seeds(k, w, num_seeds)
    calculate_entropy_vect = np.vectorize(calculate_entropy, excluded=['s_size'])
    entropies = calculate_entropy_vect(seeds, 3)
    data = pd.DataFrame(entropies, index=seeds, columns=["entropies"])
    print(seeds)
    for genome in [ecoli_k12]:  # , ecoli_BW25113, yeast]:
        hits = np.zeros(shape=(len(kmers), num_seeds), dtype=bool)
        spaced_kmers_array = determine_spaced_kmers(genome, seeds)
        for i in range(len(kmers)):
            for j in range(len(spaced_kmers_array)):
                hits[i, j] = (kmers[i] in spaced_kmers_array[j])
        print(hits)


if __name__ == "__main__":
    main()
