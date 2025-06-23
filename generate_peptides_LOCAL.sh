#!/bin/bash

# exit on error
set -e

# path to container
home_path="$HOME/Documents/PhD_Thesis/somatic_mutation_generation"
container="$home_path/peptide_netMHCpan.sif"

# input files and parameters
gtf_file="$home_path/gencode.v47.CDS.gtf"
fasta_file="$home_path/GRCh38.p14.genome.fa"
chrom="chr1"
kmer=2
pep_len=9
outdir="$home_path/results"


# path to script
script="$home_path/generate_peptides_final_version_CLUSTER.py"

# run inside container
singularity exec "$container" python "$script" \
    --gtf "$gtf_file" \
    --fasta "$fasta_file" \
    --chrom "$chrom" \
    --kmer "$kmer" \
    --pep-len "$pep_len" \
    --outdir "$outdir" \
    > "$home_path/logs/run_stdout_$chrom.log" \
    2> "$home_path/logs/run_stderr_$chrom.log"

