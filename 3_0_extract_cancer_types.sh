#!/bin/bash

# Extract cancer types and store in a Bash array
cancer_types=($(ls /home/maria/peptides/data/cbase_data_prep_MC3/output_data_preparation_MC3_original_*.txt \
  | sed -E 's|.*/output_data_preparation_MC3_original_(.*)\.txt|\1|' \
  | sort -u))

# Print the detected cancer types
echo "Detected cancer types: ${cancer_types[@]}"

# Run in parallel (change -j to control number of parallel jobs)
parallel -j 10 python3 /home/maria/peptides/scripts_human/3_a_extract_cbase_parameters.py {} ::: "${cancer_types[@]}"
