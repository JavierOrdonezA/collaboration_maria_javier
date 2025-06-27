#keeping all exons on the POSITIVE STRAND, extracting info about the start of codon poistion also
import pandas as pd
import gzip
import numpy as np
from Bio import SeqIO



#inputs
exon_coord_file = '/home/maria/peptides/auxillary/gencode_filt.bed'
fasta_file = "/home/maria/run_simulations/data/data_old/hg38.fa.gz"

#output - df with HCLCA and human exons 
output_file = '/home/maria/peptides/auxillary/exons.bed'


bases = ['A', 'C', 'G', 'T']
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def reverse_complement(seq):
    bases = list(seq)
    comp_bases=[]
    for base in bases:
        if base in complement.keys():
            comp_base = complement[base]
        else:
            comp_base = base
        comp_bases.append(comp_base)
    comp_seq = ''.join(comp_bases)
    return comp_seq[::-1]


#not using, chnaged to apply
def extract_sequences(fasta_file, coord_df, output_path):
    # Initialize a global count dictionary with only valid trinucleotides
    
    
    # Dictionary to store sequences from the FASTA file
    sequences = {}
    with gzip.open(fasta_file, "rt") as handle:
        # Parse the FASTA file and load sequences into memory
        for record in SeqIO.parse(handle, "fasta"):
            sequences[record.id] = str(record.seq)

        # List to store extracted rows for the table
    
        extracted_seqs = []
        # Extract sequences for the given regions
        for _, row in coord_df.iterrows():
            chrom = row['chr']
            start = int(row['start'])
            end = int(row['end'])
            strand = row['strand']
            

            if chrom in sequences:
                #adding flanks
                seq = sequences[chrom][start-1:end+1]
                if strand == '-':
                    seq = reverse_complement(seq)
                extracted_seqs.append(seq)
            
                    
                
            else:
                print(f"Chromosome {chrom} not found in the FASTA file.")
        coord_df['seq'] =extracted_seqs
        coord_df.to_csv(output_path, index=False)

    return 

def process_row(chr, start, end, strand, sequences):
    if chr in sequences:
        #adding flanks
        seq = sequences[chr][start-1:end+1]
        if strand == '-':
            seq = reverse_complement(seq)
    return seq


def extract_sequences_apply(fasta_file, coord_df, output_path):
    # Initialize a global count dictionary with only valid trinucleotides
    

    # Dictionary to store sequences from the FASTA file
    sequences = {}
    with gzip.open(fasta_file, "rt") as handle:
        # Parse the FASTA file and load sequences into memory
        for record in SeqIO.parse(handle, "fasta"):
            sequences[record.id] = str(record.seq)

        coord_df['seq'] = coord_df.apply(
        lambda row: process_row(row['chr'], row['start'], row['end'], row['strand'], sequences), 
        axis=1
        )

        coord_df.to_csv(output_path)

    return 



# Load the DataFrames
exon_coord_df = pd.read_csv(exon_coord_file, sep='\t', header=None, 
                            names=['chr', 'start', 'end', 'gene', 'strand'])




extract_sequences_apply(fasta_file, exon_coord_df, output_file)