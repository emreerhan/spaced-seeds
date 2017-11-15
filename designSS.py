#!/usr/bin/python
import sys
import random
import argparse
#  designSS.py
#  SSDesign
#
#  Created by Birol on 2016-03-14.
#  Modified by Emre Erhan
#  Copyright (c) 2016 Birol. All rights reserved.

# usage:
# designSS.py -n number_of_seeds -m allowed_number_of_misses -l seed_length
# rules
# m < n
# l > n


def parse_args():
    parser = argparse.ArgumentParser(
        description="Designs n spaced seeds of length l. Each bit location would have n-m seeds\
                    requiring a match. The set of seeds will be symmetric to accommodate reverse\
                    complements")
    parser.add_argument("-n", "--nseeds", type=int, default=4, help="Number of seeds (4)")
    parser.add_argument("-m", "--misses", type=int, default=2,
                        help="Number of allowed misses. Should be greater than nseeds (2)")
    parser.add_argument("-l", "--seedlen", type=int, default=100, help="Seed length. Should be more than nseeds (100)")
    parser.add_argument("-e", "--minlen", type=int, default=0, help="Minimum seed length (0)")
    parser.add_argument("-d", "--debug", action='store_true')
    args = parser.parse_args()
    return args


def design_seed(nseeds, misses, seedlen, minlen):
    current_min_length = -1
    while minlen >= current_min_length:
        # initialize seeds
        seed = [[1]*((seedlen+1)//2)+[0]*(seedlen//2) for _ in range(nseeds)]
        seed_list = range(0, nseeds)
        # scan through the seed length/2
        for i in range(0, ((seedlen+1)//2)):
            # select seeds to introduce don't care positions
            dont_care = random.sample(seed_list, nseeds-misses)
            # insert don't care positions
            for j in dont_care:
                seed[j][i] = 0
        # make seed set symmetric
        for j in range(nseeds):
            for i in range(seedlen//2):
                seed[nseeds-j-1][seedlen-i-1] = seed[j][i]

        for j in range(nseeds):
            weight = 0
            for i in range(seedlen):
                if seed[j][i] != 0:
                    weight += 1
            if weight > current_min_length:
                current_min_length = weight
        # output seeds
        for j in range(nseeds):
            weight = 0
            for i in range(seedlen):
                sys.stdout.write(str(seed[j][i]))
                if seed[j][i] != 0:
                    weight += 1
            sys.stdout.write("\t" + str(weight) + "\t" + str(seedlen-weight))
            if weight > current_min_length:
                current_min_length = weight
            print()

        for j in range(nseeds):
            for i in range(seedlen):
                sys.stdout.write(str(seed[j][i]))
            sys.stdout.write(" ")
    print()


def main():
    # get run time parameters
    args = parse_args()
    nseeds = args.nseeds  # number of seeds
    misses = args.misses  # number of allowed number of misses
    seedlen = args.seedlen  # seed length
    minlen = args.minlen
    design_seed(nseeds, misses, seedlen, minlen)


if __name__ == "__main__":
    main()
