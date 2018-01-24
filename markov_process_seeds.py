import make_seeds
import numpy as np
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Makes low entropy seeds using a Markov process. It's a dumb method...")
    # parser.add_argument("-g", "--genomes", type=str, help="Paths to the genomes, comma delimited", required=True)
    parser.add_argument("-n", "--num-seeds", type=int, help="Number of spaced seeds to generate", required=True)
    parser.add_argument("-k", "--seed-length", type=int, help="Length of seeds", required=True)
    parser.add_argument("-w", "--weight", type=int, help="Number of weighted elements per seed", required=True)
    # parser.add_argument("-t", "--transition-probability", type=float, help="Transition probability for Markov process",
    #                     required=True)
    # parser.add_argument("-e", "--entropy-bits", type=int, required=True,
    #                     help="Specify number of bits per seed for entropy calculation")
    parser.add_argument("-s", "--random-seed", type=int, help="Define a seed (default: 42)", default=42)
    parser.add_argument("-o", "--output", type=str, help="Name of output .tsv file", required=True)
    args = parser.parse_args()
    return args


def byte_to_string(byte_obj):
    return byte_obj.decode('UTF-8')


def make_entropy_seed(length, prob_transition):
    seed = np.empty(length, np.string_)
    prior = [0.5, 0.5]
    transition_prob = [[1 - prob_transition, prob_transition],
                       [prob_transition, 1 - prob_transition]]
    seed[0] = np.random.choice(['1', '0'], p=prior, size=1)[0]
    for i in range(1, length):
        last_char = int(seed[i-1])
        seed[i] = np.random.choice(['1', '0'], p=transition_prob[last_char], size=1)[0]
    byte_to_string_vect = np.vectorize(byte_to_string)
    return "".join(byte_to_string_vect(seed))


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
    temperature = 1
    while len(seeds) < num_seeds:
        prob_transition = probability_range[i]
        i += temperature
        seed = make_entropy_seed(k-2, prob_transition)
        seed = "{}{}{}".format('1', seed, '1')
        if seed.count('1') == w and seed not in seeds:
            if len(seeds) % 50 == 0:
                print('Seed: ', len(seeds))
            seeds.append(seed)
            temperature = int(1.5*temperature+100)
        if temperature > 1:
            temperature -= 1
        if i == 2000:
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
