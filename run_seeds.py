import decimal
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
        entropy -= float(count)/float(num_substrings) * math.log(decimal.Decimal(float(count)/float(num_substrings)), 2)
    return entropy


def main():

    # seed = designSS.design_seed()
    bases = ['A', 'C', 'G', 'T']
    ecoli_k12 = pyfaidx.Fasta('e_coli.fa.gz')[0][0:]
    ecoli_BW25113 = pyfaidx.Fasta('e_ecoli_BW25113.fa')[0][0:]
    yeast = pyfaidx.Fasta('Pichia_sorbitophila.fa')[0][0:]
    # print(search_sequence("AAACAAAAGTAACG", "1000011", "CGT"))
    k = 10
    w = 4
    num_kmers = 10
    kmers = []
    seeds = np.array(["".join(seed) for seed in itertools.permutations("1"*w+"0"*(k-w), k)], dtype=str)
    calculate_entropy_vect = np.vectorize(calculate_entropy, excluded=['s_size'])
    entropies = calculate_entropy_vect(seeds, 3)
    print(seeds.shape)
    print(entropies.shape)
    pd.DataFrame(entropies, index=seeds, columns=["entropies"]).to_csv("test.tsv", sep='\t')

    # for i in range(num_kmers):
    #     kmers.append(np.random.choice(bases, w, replace=True))
    #     calculate_entropy(seed, 3)


if __name__ == "__main__":
    main()
