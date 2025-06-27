from gtfparse import read_gtf
#edit to filter gecode basic gtf annotation and extract ensembl ids from my list

#rerun with full transcript, rather than basic
gtf_file = '/home/maria/data/gencode.v47.annotation.gtf.gz'
#gtf_file = '/home/maria/peptides/data/gencode.v47.basic.annotation.gtf'

gene_id_list ='/home/maria/peptides/auxillary/gtex_genes_stripped.txt'
output_path = '/home/maria/peptides/auxillary/gencode_filt.bed'

#list of expressed genes
gene_list =[]
with open(gene_id_list) as f:
    for line in f:
        gene_list.append(line.strip())

# Load the GTF file
df = read_gtf(gtf_file)
#test

print(df.columns.tolist())

#only consider protein coding genes
df_protein = df[df['gene_type']=='protein_coding']

#removing gene id version numbers, to match with expressed gene list
#filtered from 17252 pc genes to 16917
df_protein['gene_id_stripped'] = df_protein['gene_id'].str.split('.').str[0]
df_protein_filt = df_protein[df_protein['gene_id_stripped'].isin(gene_list)]

#filer quality of transcript
filt_df = df_protein_filt[df_protein_filt['transcript_support_level'].isin(['1','2'])]

# Filter for entries that have a CCDS ID
#CCDS stands for Consensus Coding Sequence — it's a high-confidence coding region annotation 
# that has been agreed upon by Ensembl, NCBI, and other reference projects
filt_df = filt_df[~filt_df['ccdsid'].isna()]

#extract only codding regions , cds plus stop codon
cds_df = filt_df[filt_df['feature'].isin(['CDS','stop_codon'])]


#ensure transcript length is a multiple of 3 for each gene
gene_groups = cds_df.groupby('gene_id')
#within each gene find cds length for each possible transcipt
transcript_gene_ids ={}
for gene_id, gene_df in gene_groups:
    #stop codon in 3utr, so add length 3 for this, plus gtf is ? based so add 1 to length
    #cds = coding seq (so only consider coding exons)
    transcript_groups = gene_df.groupby('transcript_id')
    for transcript_id, transcript_df in transcript_groups:
        #1 based so add 1 to get length of each region
        cds_length = transcript_df['end'] - transcript_df['start'] +1
        #total coding length (includes stop codon in 3utr plus cds regions)
        gene_length = cds_length.sum()
    #only valid gene annotation of multiple of 3
        if gene_length % 3 ==0:
            #using a dict like this means will only ever have one transcript selected per gene
            transcript_gene_ids[gene_id] = transcript_id



#save exon annotations for these transcript and gene ids
# transcript_gene_ids is a dict: {gene_id: transcript_id}

# Create a boolean mask using .apply()
mask = cds_df.apply(
    lambda row: transcript_gene_ids.get(row['gene_id']) == row['transcript_id'],
    axis=1
)

# Apply the mask
transcript_df = cds_df[mask]

#transcript_df = df[df['gene_id']==key and df['transcript_id']==value for key, value in transcript_gene_ids.]
#filter also based on the transcript


cut_df_exon = transcript_df[['seqname','start', 'end', 'gene_name', 'strand']]
cut_df_exon['start'] -= 1
cut_df_exon.set_index('seqname', inplace=True)
cut_df_exon.to_csv(output_path, header=False, sep ='\t')
