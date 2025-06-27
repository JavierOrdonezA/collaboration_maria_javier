INPUT="/home/maria/peptides/data/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz"

#!/bin/bash

# Input GCT file

# Output file
OUTPUT="/home/maria/peptides/auxillary/gtex_genes_filt.txt"


# Process, keeping only ensembl gene ids for genes with expression greater than 1 in atleast one sample

zcat "$INPUT" | \
    awk 'NR>2 {
        max = -1;
        for (i = 3; i <= NF; i++) {
            if ($i > max) max = $i;
        }
        if (max > 1) print $1;
    }' > "$OUTPUT"

#strip ensembl id at end
sed 's/\.[0-9]\+$//' /home/maria/peptides/auxillary/gtex_genes_filt.txt > gtex_genes_stripped.txt
