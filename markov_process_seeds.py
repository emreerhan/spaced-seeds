import make_seeds
import numpy as np
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Makes seeds using a Markov process with a dynamic transition probability to generate low to high entropy seeds.")
    # parser.add_argument("-g", "--genomes", type=str, help="Paths to the genomes, comma delimited", required=True)
    parser.add_argument("-n", "--num-seeds", type=int, help="Number of spaced seeds to generate", required=True)
    parser.add_argument("-k", "--seed-length", type=int, help="Length of seeds", required=True)
    parser.add_argument("-w", "--weight", type=int, help="Number of weighted elements per seed", required=True)
    parser.add_argument("-p", "--palindromic", help="Whether seeds should be palindromic", action='store_true', dest='palindromic')
    # parser.add_argument("-t", "--transition-probability", type=float, help="Transition probability for Markov process",
    #                     required=True)
    # parser.add_argument("-e", "--entropy-bits", type=int, required=True,
    #                     help="Specify number of bits per seed for entropy calculation")
    parser.add_argument("-s", "--random-seed", type=int, help="Define a seed (default: 42)", default=42)
    parser.add_argument("-o", "--output", type=str, help="Name of output .tsv file", required=True)
    args = parser.parse_args()
    return args


def make_entropy_seed(length, prob_transition, palindromic=False):
    even = length % 2 == 0
    if palindromic:
        length = np.floor(length/2).astype(int)
    seed = np.empty(length, 'S1')
    prior = [0.5, 0.5]
    transition_prob = [[1 - prob_transition, prob_transition],
                       [prob_transition, 1 - prob_transition]]
    seed[0] = np.random.choice(['1', '0'], p=prior, size=1)[0]
    for i in range(1, length):
        last_char = int(seed[i-1])
        seed[i] = np.random.choice(['1', '0'], p=transition_prob[last_char], size=1)[0]
    if palindromic:
        seed_reverse = np.flipud(seed)
        if not even:
            seed = np.append(seed, [0])
        seed = np.append(seed, seed_reverse)
    return seed.tostring().decode('utf-8')


def main():
    args = parse_args()
    k = args.seed_length
    w = args.weight
    # prob_transition = 1-args.transition_probability
    num_seeds = args.num_seeds
    np.random.seed(args.random_seed)
    prob_min = 0.001
    prob_max = 0.999
    probability_range = np.linspace(prob_min, prob_max, 2000)
    seeds = [] 
    print('Generating {} seeds with k = {}, w = {}, p(transition) ∈ [{}, {}]'.format(
        num_seeds, k, w, prob_min, prob_max))
    i = 0
    # temperature = 1
    while len(seeds) < num_seeds:
        prob_transition = probability_range[i]
    #    i += temperature
        seed = make_entropy_seed(k-2, prob_transition, palindromic=args.palindromic)
        seed = "{}{}{}".format('1', seed, '1')
        if seed.count('1') == w and seed not in seeds:
            if len(seeds) % 50 == 0:
                print('Seed: ', len(seeds))
            seeds.append(seed)
        if i == len(probability_range):
            i = 0
    calculate_entropy_vect = np.vectorize(make_seeds.calculate_entropy, excluded=['s_size'])
    entropies_1 = calculate_entropy_vect(seeds, 2)
    entropies_2 = calculate_entropy_vect(seeds, 3)
    entropies_3 = calculate_entropy_vect(seeds, 4)
    data = pd.DataFrame({'2bit': entropies_1, '3bit': entropies_2, '4bit': entropies_3},
                        index=seeds,
                        columns=['2bit', '3bit', '4bit'])
    data.to_csv(args.output, sep='\t')
    print(data.describe())


if __name__ == "__main__":
    main()
