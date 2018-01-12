import numpy as np
import pandas as pd
import math
import pyfaidx
import argparse
import pickle
# import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Runs a number of spaced seeds against specified genomes
                    to determine if seed entropy affects specificity.""")
    parser.add_argument("-g", "--genomes", type=str, help="Paths to the genomes, comma delimited", required=True)
    parser.add_argument("-n", "--num-seeds", type=int, help="Number of spaced seeds to sample", required=True)
    parser.add_argument("-s", "--seeds", type=str, help="Seeds tsv file, output of make_seeds.py", required=True)
    # parser.add_argument("-o", "--output", type=str, help="Name of output .tsv file", required=True)
    args = parser.parse_args()
    return args


def sample_seeds(seeds, entropies, sample_size):
    select_indexes = np.logspace(0, math.log10(len(entropies)-1), sample_size).astype(int)
    sorted_indexes = np.argsort(entropies)
    indexes = np.sort(sorted_indexes[select_indexes])
    return np.array(seeds)[indexes]


def determine_word_frequencies(genome_sequence, spaced_seed):
    k = len(spaced_seed)
    weighted_indexes = []
    kmers = {}
    for i in range(k):
        if spaced_seed[i] == '1':
            weighted_indexes.append(i)
    w = len(weighted_indexes)
    for i in range(len(genome_sequence) - k + 1):
        word = ['']*k
        for j, weighted_index in zip(range(w), weighted_indexes):
            word[j] = genome_sequence[i:i+k][weighted_index]
        word_str = "".join(word)
        kmers[word_str] = kmers.get(word_str, 0) + 1
    return kmers


def determine_word_uniqueness(words_dict):
    np_values = np.array(list(words_dict.values()))
    unique_items_index = np_values == 1
    p_uniqueness = len(np_values[unique_items_index]) / len(np_values)
    return p_uniqueness


def main():
    args = parse_args()
    genome = str(pyfaidx.Fasta(args.genomes)[0][0:])
    num_seeds = args.num_seeds
    data = pd.read_csv(args.seeds, sep='\t', index_col=0)
    seeds = data.index.values
    prefix = "e_coli"
    seed_sample = sample_seeds(seeds, data['3bit'].values, num_seeds)
    determine_words_vect = np.vectorize(determine_word_frequencies, excluded=['genome_sequence'])
    words_list = determine_words_vect(genome, seed_sample)
    deter_uniqueness_vect = np.vectorize(determine_word_uniqueness)
    uniquenesses = deter_uniqueness_vect(words_list)
    print("seed", "uniqueness")
    for seed, uniqueness in zip(seed_sample, uniquenesses):
        print(seed, uniqueness, sep='\t')

    # for idx, words in enumerate(words_list):
    #     output = open('{}/seed{}.pkl'.format(prefix, idx), 'wb')
    #     pickle.dump(words, output)
    #     output.close()
    # seed_output = open('{}/seeds.txt'.format(prefix), 'w')
    # for idx, seed in enumerate(seed_sample):
    #     print(idx, '\t', seed, file=seed_output)
    # seed_output.close()

    # out_data = pd.DataFrame.from_dict(words_list.tolist())
    # out_data.index = seed_sample
    # out_data.T.to_csv(args.output, sep='\t')


if __name__ == "__main__":
    main()
