# Spaced Seed Design
There are many heuristics for spaced seed design. However, most of these methods are not good for designing large spaced seeds (around k=60 and w=20).

My undergraduate directed studies project examined if Shannon entropy can be used as an approximation for spaced seed quality for database searching, particularily when using [BioBloom tools](https://github.com/bcgsc/biobloom).

## Scripts
* make_seeds.py: Randomly generates spaced seeds of a given k and w
* markov_process_seeds.py: Generate spaced seeds of varying entropy for a given k and w
* determine_uniqueness.py: Determines the uniqueness of the set of words produced by a tsv of spaced seeds for a given genome.
* select_multi_spaced_seeds.py: Generates a list of multiple spaced seeds, where each set has 5 spaced seeds, designed for use in BBT.

## Manuscript
https://goo.gl/Qaed8m
