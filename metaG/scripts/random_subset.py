#!/usr/bin/env python

import random
random.seed(42)
random_rows = random.sample(range(snakemake.params.total), k=snakemake.params.number) # change if needed
with open(snakemake.input[0], "r") as ifile, open(snakemake.output[0], "w") as ofile:
    count = 0
    for line in ifile:
        if count == 0 or (count-1) in random_rows:
            ofile.write(line)
        count += 1
