import numpy as np
import pandas as pd
import argparse
from itertools import combinations


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate sets of good, average, and bad multi-spaced")
    parser.add_argument("-s", "--seeds", type=str, help="Seeds tsv file, output of make_seeds.py", required=True)
    parser.add_argument("-p", "--palindromes", type=str, help="Palindromic seeds tsv file, output of make_seeds.py",
                        required=True)
    # parser.add_argument("-o", "--output", type=str, help="Name of output .tsv file", required=True)
    args = parser.parse_args()
    return args


def overlap_complexity(seeds):
    # For now assume 2 seeds
    complexity = 0
    for seed1, seed2 in combinations(seeds, 2):
        for i in range(1, len(seed1)):
            sigma1 = get_sigma(seed1[:i], seed2[len(seed2)-i:])
            sigma2 = get_sigma(seed1[len(seed1)-i:], seed2[:i])
            complexity += 2 ** sigma1 + 2 ** sigma2
        final_sigma = get_sigma(seed1, seed2)
        complexity += 2 ** final_sigma
    return complexity


def get_sigma(s1, s2):
    sigma = 0
    for i, j in zip(s1, s2):
        if i == '1' and j == '1':
            sigma += 1
    return sigma


def main():
    args = parse_args()
    seed_data = pd.read_csv(args.seeds, sep='\t', index_col=0)
    palindrome_seeds = pd.read_csv(args.palindromes, sep='\t', index_col=0)
    # seeds = seed_data.index.values
    num_seeds = len(seed_data.index)
    seed_data = seed_data.sort_values(by='2bit', ascending=False)
    palindrome_data = palindrome_seeds.sort_values(by='2bit', ascending=False)
    seed_index = np.linspace(0, len(seed_data.index.values)-1, 100).astype(int)
    seed_combos = combinations(seed_data.index.values[seed_index], 2)
    seed_combos = np.array([(*combo, combo[0][::-1], combo[1][::-1]) for combo in seed_combos])
    new_seed_combos = []
    for seed_combo in seed_combos:
        for palindrome in palindrome_data.index.values[:5]:
            new_seed_combos.append([*seed_combo, palindrome])

    overlap_complexities = np.empty(len(new_seed_combos))
    for i, combo in enumerate(new_seed_combos):
        overlap_complexities[i] = overlap_complexity(combo)
    for seeds, complexity in zip(new_seed_combos, overlap_complexities):
        print(*seeds, sep=' ', end='\t')
        print(complexity)


if __name__ == '__main__':
    main()
