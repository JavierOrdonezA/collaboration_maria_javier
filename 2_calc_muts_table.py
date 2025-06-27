#finding subs at all poistions plus recording the sequence positoins when merging exons 
import pandas as pd
from pathlib import Path



#usage, input is exon_df, contains flanking introns, and duplicate regions subtracted
exon_file = '/home/maria/peptides/auxillary/exons.bed'
gene_annotations = '/home/maria/peptides/auxillary/gencode_filt.bed'


#output
output = '/home/maria/peptides/auxillary/possible_muts'



#functions required for finding subs

bases = ['A', 'C', 'G', 'T']
trinucleotides = [a + b + c for a in bases for b in bases for c in bases]



def add_trinucs_remove_flanks(reconstructed_gene_df, trinucleotides):
    reconstructed_gene_df['trinucs'] = reconstructed_gene_df['seq'].apply(
        lambda x: list_trinucs_of_flanked_seq(x, trinucleotides)
    )
    #REMOVE FLANKS
    reconstructed_gene_df['seq'] = reconstructed_gene_df['seq'].apply(
        lambda x: x[1:-1]
    )
    return reconstructed_gene_df


#fills gaps in exon sequences in df, checked when start, end is a subinterval of exon start, end (by intersection should always be true)
def fill_gaps_in_exons(exon_combined_df):
    # Calculate how many Ks to add at the beginning
    pad_start = (exon_combined_df['start'] - exon_combined_df['exon_start'])
    
    # Calculate how many Ks to add at the end
    pad_end = (exon_combined_df['exon_end'] - exon_combined_df['end'])

    # Generate padded sequences
    k_start = pad_start.apply(lambda x: 'K' * x)
    k_end = pad_end.apply(lambda x: 'K' * x)
    k_start_trinuc = pad_start.apply(lambda x: [-1] * x)
    k_end_trinuc = pad_end.apply(lambda x: [-1] * x)

    # Update the sequence column with padded values
    exon_combined_df['seq'] = k_start + exon_combined_df['seq'] + k_end
    # Update the trinucs column with padded values
    exon_combined_df['trinucs'] = k_start_trinuc + exon_combined_df['trinucs'] + k_end_trinuc
  
    exon_combined_df.drop(columns=['start', 'end'], axis =1, inplace=True)
    return exon_combined_df



def add_extra_exons(exon_combined_df, exon_coordinates_df):
    
    # Step 1: Find missing exons — those in exon_coordinates_df not already in exon_combined_df
    merged = exon_coordinates_df.merge(
        exon_combined_df[['chr', 'exon_start', 'exon_end','gene','strand']],
        on=['chr', 'exon_start', 'exon_end','gene','strand'],
        how='left',
        indicator=True
    )

    #print(merged)
    # Step 2: Filter only missing ones
    missing_exons_df = merged[merged['_merge'] == 'left_only'].copy()
    # Step 3: Add dummy sequence of 'K's of appropriate length
    missing_exons_df['seq'] = missing_exons_df.apply(
        lambda row: 'K' * (row['exon_end'] - row['exon_start']), axis=1
    )
    missing_exons_df['trinucs'] = missing_exons_df.apply(
        lambda row: [-1] * (row['exon_end'] - row['exon_start']), axis=1
    )
    
    missing_exons_df.drop(columns=['_merge'],inplace=True)
    # Step 4: Combine the original and missing exons
    full_df = pd.concat([exon_combined_df, missing_exons_df], ignore_index=True)

    return full_df

#assuming zero based coordinates
def list_trinucs_of_flanked_seq(seq,trinucleotides):
    trinucs = []
    
    k = len(seq)
    for i in range(1, k-1):
      #iterating through all positions, exluding flanking bases 
        trinuc_string = seq[i-1:i+2]
        if trinuc_string in trinucleotides:
            trinuc = trinucleotides.index(trinuc_string)
            trinucs.append(trinuc)  
        else:
            trinucs.append(-1)
        
    return trinucs


def add_seq_pos(exon_start,exon_end):
    return list(range(exon_start, exon_end))



from itertools import chain

#editted here assuming REFERENCE HUMAN, ie do not need to add gaps etc, start =exon start
#also extract a list of seq positions here 
def extract_gene_seq_trinucs(gene_df):
    gene_df['pos'] = gene_df.apply(
        lambda row: add_seq_pos(row['start'],row['end']), axis=1
    )
    gene_df.to_csv('/home/maria/target_size_MRCA/test')
    strand = gene_df.iloc[0]['strand']
    if strand == '+':
        gene_df = gene_df.sort_values(by='start').reset_index(drop=True)
    else:
        gene_df = gene_df.sort_values(by='end', ascending=False).reset_index(drop=True)
    seq = ''.join(gene_df['seq'].tolist())
    trinucs = list(chain.from_iterable(gene_df['trinucs']))
    
    pos = list(chain.from_iterable(gene_df['pos']))
    
    return seq, trinucs, pos




#define the codon-to-amino acid mapping
codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}

#function to check what type of mutation has occured
def find_mutation_type(original_codon, mutated_codon):
    # Translate both codons using the codon table
    if (original_codon in codon_table) and (mutated_codon in codon_table):
        original_aa = codon_table.get(original_codon)
        mutated_aa = codon_table.get(mutated_codon)
        if original_aa == mutated_aa:
            type = 0  #mutatoin is syn
        else:
            type = 1 #mutation is non-syn
    else:
        type = "not a valid codon"
    return type


def generate_mutations(seq, trinucs, trinucleotides, seq_positions):

    k = len(seq)
    context_indexes =[]
    refs =[]
    alts =[]
    mutation_type = []
    indexes =[]
    #finding start index of first codon
    i=0
    while i<(k-2): #up to last codon
      
        ref_codon = seq[i:i+3]
        
        if ref_codon in trinucleotides:
            for m in range(3):
                seq_pos = i+m
                og_nuc = seq[seq_pos]
                possible_muts = ['A', 'C', 'G', 'T']
                possible_muts.remove(og_nuc) 
                
                
                    #context_index = find_context_index(ref_trinucleotide)
                for alt in possible_muts:
                # find mutation type
                    if m == 0:
                        alt_codon =  alt + ref_codon[1] + ref_codon[2]
                    if m == 1:
                        alt_codon = ref_codon[0] + alt + ref_codon[2]
                    if m == 2:
                        alt_codon = ref_codon[0] + ref_codon[1] + alt

                    type = find_mutation_type(ref_codon, alt_codon)
                
                    
                        
                    ref_trinucleotide = trinucs[seq_pos]

                    
                    if ref_trinucleotide != -1:
                        alts.append(alt)
                        context_indexes.append(ref_trinucleotide)
                        mutation_type.append(type)
                        refs.append(og_nuc)
                        indexes.append
                        indexes.append(seq_positions[seq_pos])
                        
        i=i+3
    return context_indexes, refs, alts, mutation_type, indexes


reconstructed_df = pd.read_csv(exon_file)
coord_df = pd.read_csv(gene_annotations, delimiter='\t', names = ['chr','exon_start', 'exon_end','gene','strand'])


reconstructed_trinuc_df = add_trinucs_remove_flanks(reconstructed_df, trinucleotides)

'''
full_exons_df = fill_gaps_in_exons(reconstructed_trinuc_df)
print(full_exons_df)
all_exons_df = add_extra_exons(full_exons_df, coord_df)
print(all_exons_df)
'''

all_exons_df = reconstructed_trinuc_df
#correct
f = open(output, "w")

for gene in coord_df['gene'].unique():

    gene_df = all_exons_df[all_exons_df['gene']==gene]
    
    seq,trinucs, pos =extract_gene_seq_trinucs(gene_df)
    if len(seq) %3==0:
        chr = gene_df['chr'].iloc[0]
        
        context_indexes, refs, alts, mutation_type, indexes = generate_mutations(seq, trinucs, trinucleotides, pos)
        for k in range(len(alts)):
            f.write(chr + "\t" + str(indexes[k]) + "\t" + gene + "\t" + str(mutation_type[k]) + "\t" + refs[k] + "\t" + alts[k] + "\t" + str(context_indexes[k])+ "\n")


f.close()
        
       
 


















        
                
                
                
                






                

        


