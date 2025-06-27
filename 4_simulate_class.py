import numpy as np
import pandas as pd
import sys
from NeutralEvolutionSimulator import NeutralEvolutionSimulator

cancer_type = sys.argv[1]

simulator = NeutralEvolutionSimulator(
    n_subs=100000,
    coord_file='/home/maria/peptides/auxillary/exons.bed',
    cbase_output=f'/home/maria/peptides/auxillary/cbase/cbase_filt_{cancer_type}.csv',
    syn_probs='/home/maria/peptides/auxillary/syn_probs.csv',
    non_syn_probs='/home/maria/peptides/auxillary/non_syn_probs.csv',
    output=f'/home/maria/peptides/output/simulated_exons_flank_{cancer_type}'
)

simulator.run()
