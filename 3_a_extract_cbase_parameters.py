import pandas as pd
import sys

cancer_type = sys.argv[1]

data_preparation_file=f"/home/maria/peptides/data/cbase_data_prep_MC3/output_data_preparation_MC3_original_{cancer_type}.txt"
qvalues_file=f"/home/maria/target_size_mammals_gene_orthogs/target_size_human_longitudinal/data_cbase/q_values_MC3_original_{cancer_type}.txt"
# Output file
output_file=f"/home/maria/peptides/auxillary/cbase/cbase_filt_{cancer_type}.csv"

qvalues_df = pd.read_csv(qvalues_file, skiprows = 1, index_col =0, delimiter='\t')

data_df = pd.read_csv(data_preparation_file, index_col =0, delimiter='\t')

merged_df = qvalues_df[['d(m+k)/ds']].join(data_df)
merged_df['lambda_n'] = merged_df['lambda_s'] * (merged_df['l_k'] + merged_df['l_m']) /merged_df['l_s']
merged_df['dN_lambda_n'] = merged_df['lambda_n'] * merged_df['d(m+k)/ds']

cut_df = merged_df[['dN_lambda_n','lambda_s']]
cut_df.to_csv(output_file)