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
    # select_indexes = np.logspace(0, math.log10(len(entropies)-1), sample_size).astype(int)
    select_indexes = np.linspace(0, len(entropies)-1, sample_size, dtype=int) #.astype(int)
    sorted_indexes = np.argsort(entropies)
    indexes = np.sort(sorted_indexes[select_indexes])
    return np.array(seeds)[indexes]


def reverse_complement(numpy_sequence):
    complement = {b'A': b'T',
                  b'T': b'A',
                  b'C': b'G',
                  b'G': b'C'}
    rev_seq = np.flipud(numpy_sequence)
    return np.vectorize(lambda nt: complement[nt])(rev_seq)


def determine_word_frequencies(genome_sequence, spaced_seed):
    spaced_seed = str(spaced_seed)
    k = len(spaced_seed)
    weighted_indexes = []
    kmers = {}
    for i in range(k):
        if spaced_seed[i] == '1':
            weighted_indexes.append(i)
    w = len(weighted_indexes)
    for i in range(len(genome_sequence) - k + 1):
        word = np.empty(w, 'S1')
        for j, weighted_index in zip(range(w), weighted_indexes):
            word[j] = genome_sequence[i:i+k][weighted_index]
        word_str = word.tostring() #.decode('utf-8')
        reverse_word_str = reverse_complement(word).tostring() #.decode('utf-8')
        canonical_word = min(word_str, reverse_word_str)
        # kmers[canonical_word] = False if canonical_word in kmers else True
        kmers[canonical_word] = kmers[canonical_word] + [i] if canonical_word in kmers else [i]
    return kmers


def determine_unique_frames(kmers, genome_length, k):
    num_frames = genome_length - k + 1
    frames_count = np.empty(num_frames, np.int8)
    for frames in kmers.values():
        for frame in frames:
            if frames_count[frame] < 3:
                frames_count[frame] += 1
    return frames_count == 1


def determine_word_uniqueness(words_dict):
    np_values = np.array(list(words_dict.values()))
    # unique_items_index = np_values == True
    p_uniqueness = len(np_values[np_values]) / len(np_values)
    return p_uniqueness


def main():
    args = parse_args()
    genome = str(pyfaidx.Fasta(args.genomes)[0][0:])
    num_seeds = args.num_seeds
    data = pd.read_csv(args.seeds, sep='\t', index_col=0)
    seeds = data.index.values
    k = 60
    prefix = "e_coli"
    seed_sample = sample_seeds(seeds, data['3bit'].values, num_seeds)
    # seed_sample = seeds
    determine_words_vect = np.vectorize(determine_word_frequencies, excluded=['genome_sequence'])
    kmers_list = determine_words_vect(genome, seed_sample)
    # kmers = determine_word_frequencies(genome, seed_sample[0])
    determine_unique_frames_vect = np.vectorize(determine_unique_frames, excluded=['genome_length', 'k'])
    print(determine_unique_frames(kmers_list[0], len(genome), k))
    num_frames = len(genome) - k + 1
    # unique_frames_list = determine_unique_frames_vect(kmers_list, len(genome), k)
    unique_frames_list = np.empty((len(kmers_list), num_frames), dtype=np.bool)
    for idx, kmers in enumerate(kmers_list):
        unique_frames_list[idx] = determine_unique_frames(kmers, len(genome), k)
    first_choices = np.random.choice(len(unique_frames_list), 10)
    second_choices = np.random.choice(len(unique_frames_list), 10)
    combined_uniqueness = np.logical_or(unique_frames_list[first_choices], unique_frames_list[second_choices])
    unique_frames = sum(combined_uniqueness)

    # deter_uniqueness_vect = np.vectorize(determine_word_uniqueness)
    # uniquenesses = deter_uniqueness_vect(words_list)
    print("seed", "uniqueness")
    for seed1, seed2, uniqueness in zip(seed_sample[first_choices], seed_sample[second_choices], unique_frames):
        print(seed1, seed2, uniqueness, sep='\t')


if __name__ == "__main__":
    main()
