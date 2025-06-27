#extracting lists of trinucleotide context probabilties from SBS5 file
#possible muts file is extremly large s only saving one per gene

import pandas as pd

intron_matrix_file = '/home/maria/target_size_mammals_gene_orthogs/signatures/TSB_signatures_scaled/SBS5'
possible_subs_file = '/home/maria/peptides/auxillary/possible_muts'
output_syn = '/home/maria/peptides/auxillary/syn_probs.csv'
output_non_syn = '/home/maria/peptides/auxillary/non_syn_probs.csv'

bases = ['A', 'C', 'G', 'T']
trinucleotides = [a + b + c for a in bases for b in bases for c in bases]
intron_matrix_df = pd.read_csv(intron_matrix_file, index_col=0)

def extract_trio_prob(trinuc_string, alt):
    return intron_matrix_df.at[trinuc_string, alt]

# Write headers once
with open(output_syn, 'w') as f:
    f.write("gene,pos,prob\n")
with open(output_non_syn, 'w') as f:
    f.write("gene,pos,prob\n")

# Chunked reading and buffering
chunk_reader = pd.read_csv(
    possible_subs_file,
    sep='\t',
    header=None,
    names=['chr', 'pos', 'gene', 'mut_type', 'ref', 'alt', 'trinuc_context'],
    chunksize=100000
)

gene_buffer = pd.DataFrame()

for chunk in chunk_reader:
    # Combine current chunk with leftover rows from previous chunk
    combined = pd.concat([gene_buffer, chunk], ignore_index=True)

    # Get all genes present and sorted
    genes_in_chunk = combined['gene'].unique()

    # Assume the **last gene** may be incomplete — hold it back
    complete_genes = genes_in_chunk[:-1]
    incomplete_gene = genes_in_chunk[-1]

    # Process complete genes only
    for gene in complete_genes:
        gene_df = combined[combined['gene'] == gene]
        if gene_df.empty:
            continue

        for mut_type, out_file in [(0, output_syn), (1, output_non_syn)]:
            filtered = gene_df[gene_df['mut_type'] == mut_type]
            if filtered.empty:
                continue
            try:
                probs = filtered.apply(
                    lambda row: extract_trio_prob(
                        trinucleotides[row['trinuc_context']],
                        row['alt']
                    ),
                    axis=1
                )
                #normalising per gene probabilities
                norm_probs = probs / probs.sum()
                filtered = filtered.assign(prob=norm_probs)
                filtered[['gene', 'pos', 'alt', 'prob']].to_csv(out_file, mode='a', index=False, header=False)
            except KeyError:
                print(f"Skipped gene {gene} (type {mut_type}): trinuc or alt not found")

    # Save incomplete gene to buffer for next chunk
    gene_buffer = combined[combined['gene'] == incomplete_gene]

# Process any leftover gene in final buffer
if not gene_buffer.empty:
    gene = gene_buffer['gene'].iloc[0]
    for mut_type, out_file in [(0, output_syn), (1, output_non_syn)]:
        filtered = gene_buffer[gene_buffer['mut_type'] == mut_type]
        if filtered.empty:
            continue
        try:
            probs = filtered.apply(
                lambda row: extract_trio_prob(
                    trinucleotides[row['trinuc_context']],
                    row['alt']
                ),
                axis=1
            )
            norm_probs = probs / probs.sum()
            filtered = filtered.assign(prob=norm_probs)
            filtered[['gene', 'pos', 'alt', 'prob']].to_csv(out_file, mode='a', index=False, header=False)
        except KeyError:
            print(f"Skipped final gene {gene} (type {mut_type})")
