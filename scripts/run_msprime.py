#!/usr/bin/env python3
import msprime as msp
import sys
import pickle

with open(sys.argv[1], 'rb') as f:
    kwargs, contig_id, out = pickle.load(f)
sim = msp.simulate(**kwargs)
sim.write_vcf(open(out, "wt"), ploidy=2, contig_id=contig_id)
