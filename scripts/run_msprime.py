#!/usr/bin/env python3
import msprime
import pickle
import sys

kwargs, contig_id, orig = pickle.load(open(sys.argv[1], "rb"))
sim = msprime.simulate(**kwargs)
sim.write_vcf(open(orig, "wt"), ploidy=2, contig_id=contig_id)
