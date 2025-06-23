# %%
"""
Genomics Analysis Pipeline for Isoform Selection and Mutation Analysis

This pipeline performs:
1. Isoform selection and canonical transcript identification
2. Sequence extraction with context padding
3. K-mer analysis and mutation generation
4. Peptide sequence analysis with mutations

Author: J.Ordonez A

"""

from pathlib import Path
from typing import List, Tuple, Dict, Optional, Any

import pandas as pd
import numpy as np
from pyfaidx import Fasta
from gtfparse import read_gtf
from tqdm import tqdm


import sys
import logging
import warnings
warnings.filterwarnings("ignore")

# %%
"""
# ) Extract only CDS lines from the GTF
with open("gencode.v47.annotation.gtf") as gtf, 
        open("gencode.v47.CDS.gtf", "w") as out:
    for line in gtf:
        if line.startswith("#"):
            out.write(line)
        else:
            fields = line.split("\t")
            if fields[2] == "CDS":
                out.write(line)

with open("gencode.v47.CDS.gtf") as f:
    lines = f.read().splitlines()
data_lines = [line for line in lines if not line.startswith("#")]
records = [line.split("\t") for line in data_lines]
cds_records = [fields for fields in records if fields[2] == "CDS"]
"""


class GenomicsConfig:
    # Configuration class for genomics analysis parameters.#

    def __init__(self):
        # File paths - update these to your actual paths
        self.path_cds_annotation = '/Users/fordonez/Documents/PhD_Thesis/' \
            + 'somatic_mutation_generation/gencode.v47.CDS.gtf'
        self.path_grch38 = '/Users/fordonez/Documents/PhD_Thesis/somatic_mutation_generation/' \
            + 'GRCh38.p14.genome.fa'

        # Analysis parameters
        self.k_mer_neighbors = 2  # Padding for k-mer context
        self.peptide_length = 9   # Number of amino acids in peptide
        self.window_length = 3 * self.peptide_length  # DNA window length

        # Chromosome to analyze
        self.chromosome = 'chr1'


class CodonTable:
    """Standard genetic code codon table."""

    TABLE = {
        'TAA': '*', 'TAG': '*', 'TGA': '*',  # Stop codons
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',  # Alanine
        'TGC': 'C', 'TGT': 'C',  # Cysteine
        'GAC': 'D', 'GAT': 'D',  # Aspartic acid
        'GAA': 'E', 'GAG': 'E',  # Glutamic acid
        'TTC': 'F', 'TTT': 'F',  # Phenylalanine
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',  # Glycine
        'CAC': 'H', 'CAT': 'H',  # Histidine
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I',  # Isoleucine
        'ATG': 'M',  # Methionine (start codon)
        'AAA': 'K', 'AAG': 'K',  # Lysine
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'TTA': 'L', 'TTG': 'L',  # Leucine
        'AAC': 'N', 'AAT': 'N',  # Asparagine
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',  # Proline
        'CAA': 'Q', 'CAG': 'Q',  # Glutamine
        'AGA': 'R', 'AGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',  # Arginine
        'AGC': 'S', 'AGT': 'S', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',  # Serine
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',  # Threonine
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',  # Valine
        'TGG': 'W',  # Tryptophan
        'TAC': 'Y', 'TAT': 'Y',  # Tyrosine
    }

    @classmethod
    def translate(cls, codon: str) -> str:
        """Translate a codon to amino acid."""
        return cls.TABLE.get(codon.upper(), 'X')

    @classmethod
    def translate_sequence(cls, dna_seq: str) -> str:
        """Translate DNA sequence to protein sequence."""
        if len(dna_seq) % 3 != 0:
            warnings.warn(
                f"DNA sequence length {len(dna_seq)} is not divisible by 3")

        protein = []
        for i in range(0, len(dna_seq), 3):
            codon = dna_seq[i:i+3].upper()
            if len(codon) == 3:
                protein.append(cls.translate(codon))
        return "".join(protein)


class IsoformSelector:
    """Handles isoform selection and canonical transcript identification."""

    @staticmethod
    def get_isoforms_by_gene(df: pd.DataFrame, chromosome: str) -> pd.DataFrame:
        """
        Given a GTF‐loaded DataFrame with separate gene_id, transcript_id, gene_name columns,
        return one row per gene on `chromosome` with:
        - gene_id
        - gene_name
        - n_isoforms   (number of distinct transcript_id)
        - isoforms     (sorted list of transcript_id)
        """
        # 1) Filter to the chromosome of interest
        df_chr = df[df['seqname'] == chromosome]

        # 2) Group by gene_id + gene_name
        summary = (
            df_chr
            .groupby(['gene_id', 'gene_name'])['transcript_id']
            .agg(
                n_isoforms='nunique',
                isoforms=lambda txs: sorted(txs.unique())
            )
            .reset_index()
        )
        return summary

    @staticmethod
    def filter_isoforms_by_phase(df_annotations: pd.DataFrame, chromosome: str) -> pd.DataFrame:
        """
        Filter an annotation DataFrame for a given chromosome and keep only those
        transcript_id (within each gene_id) whose total segment_length is a multiple
        of 3, after applying:

        1) Chromosome selection.
        2) CCDS flag extraction from the 'tag' column.
        3) Quality filtering:
            - protein_coding must have CCDS=True
            - protein_coding must have MANE_Select=True
            - protein_coding must have level in {1,2}
        4) Calculation of segment_length = end - start + 1

        Parameters
        ----------
        df_annotations : pd.DataFrame
            Must contain columns:
            ['seqname', 'tag', 'gene_type', 'level', 'start', 'end',
            'gene_id', 'transcript_id'].
        chromosome : str
            Chromosome to filter on (e.g. 'chr1').

        Returns
        -------
        pd.DataFrame
            Subset of the original DataFrame containing only those rows whose
            transcript_id pass the “mod-3 = 0” test.
        """

        # 1) copy & chromosome filter
        df = df_annotations.copy()
        df_chr = df[df['seqname'] == chromosome].copy()

        # 2)Extract CCDS flag
        df_chr['CCDS'] = df_chr['tag'].str.contains(r'\bCCDS\b', na=False)

        df_chr['MANE_Select'] = df_chr['tag'].str.contains(
            r'\bMANE_Select\b', na=False)

        # 3) Apply quality filters
        is_pc = df_chr['gene_type'] == 'protein_coding'
        # pass_pc = is_pc & df_chr['CCDS']
        pass_pc = is_pc & (df_chr['CCDS'] | df_chr['MANE_Select'])
        # pass_npc = ~is_pc & df_chr['level'].isin([1, 2])
        pass_npc = is_pc & df_chr['level'].isin([1, 2])

        df_q = df_chr[pass_pc | pass_npc].copy()

        # 4) Calculate segment lengths
        df_q['segment_length'] = df_q['end'] - df_q['start'] + 1

        # 5) sum lengths per transcript, keep those %3 == 0
        seg_sum = (
            df_q
            .groupby(['gene_id', 'transcript_id'])['segment_length']
            .sum()
            .reset_index(name='total_length')
        )
        valid = seg_sum[seg_sum['total_length'] %
                        3 == 0][['gene_id', 'transcript_id']]

        # join back to get only valid rows
        df_valid = df_q.merge(
            valid, on=['gene_id', 'transcript_id'], how='inner')

        # drop the helper column
        df_valid = df_valid.drop(columns=['segment_length'])

        df_valid["exon_number"] = pd.to_numeric(
            df_valid["exon_number"], errors="coerce")
        return df_valid

    @staticmethod
    def pick_canonical_transcript_2(df: pd.DataFrame) -> pd.DataFrame:

        # From a GTF‐style DataFrame with columns including:
        # ['gene_id','transcript_id','feature','start','end','tag']
        # return one row per gene giving the chosen canonical transcript_id.

        df = df.copy()

        # 1) Calculate CDS lengths
        cds = df[df['feature'] == 'CDS'].copy()
        cds['cds_len'] = cds['end'] - cds['start'] + 1
        cds_len = cds.groupby('transcript_id')[
            'cds_len'].sum().rename('total_cds_len')

        # 2) Extract boolean flags for each priority tag
        priority_tags = ['MANE_Select', 'appris_principal_1', 'Ensembl_canonical',
                         'CCDS']

        for tag in priority_tags:
            df[tag] = df['tag'].fillna('').str.contains(
                rf'\b{tag}\b', na=False)

        # 3) Build a per-transcript summary
        tx = (
            df[['gene_id', 'transcript_id', 'gene_name']]
            .drop_duplicates()
            .set_index('transcript_id')
        )
        tx = tx.join(cds_len, how='left').fillna({'total_cds_len': 0})

        for tag in priority_tags:
            tx[tag] = df.drop_duplicates(subset=['transcript_id'])[tag].values

        # 4) Define a scoring function per transcript
        def calculate_priority_score(row):
            priorities = {
                'MANE_Select': 0,
                'CCDS': 1,
                'appris_principal_1': 2,
                'Ensembl_canonical': 3  # protein functionality score
            }
            for tag, priority in priorities.items():
                if row[tag]:
                    return priority
            return 4  # No priority tag found

        tx['priority_score'] = tx.apply(calculate_priority_score, axis=1)

        # Select best transcript per gene
        # 5) For each gene, pick the transcript with
        #   a) lowest priority_score
        #   b) if tie, highest total_cds_len
        #   c) if still tie, choose the first element (oldest transcript_id)
        def select_best_transcript(sub_df):
            # a) filter by min score
            best = sub_df[sub_df['priority_score']
                          == sub_df['priority_score'].min()]
            # b) if ties, filter by max length
            best = best[best['total_cds_len'] == best['total_cds_len'].max()]
            # c) if still ties, pick first by transcript_id sort
            return best.sort_index().iloc[0]

        canonical = (
            tx
            .reset_index()
            .groupby('gene_id')
            .apply(select_best_transcript)
            .reset_index(drop=True)
            [['gene_id', 'transcript_id', 'gene_name',
                'priority_score', 'total_cds_len']]
        )

        return canonical

    @staticmethod
    def pick_canonical_transcript(df: pd.DataFrame) -> pd.DataFrame:
        """
        From a GTF‐style DataFrame with columns
        ['gene_id', 'transcript_id', 'gene_name', 'feature', 'start', 'end', 'tag'],
        return one row per gene with the chosen canonical transcript.

        Selection steps per gene:
        1) Compute total CDS length per transcript.
        2) Extract boolean flags for each priority tag.
        3) Summarize one row per transcript.
        4) Within each gene, pick the best transcript by:
            a) Lowest priority_score (0=MANE_Select, 1=CCDS, 2=appris_principal_1,
                3=Ensembl_canonical, 4=none).
            b) Hierarchical tag filtering: for each tag in order, if it splits the
                set (some have it, some don’t), keep only those that have it.
            c) If still tied, highest total CDS length.
            d) If still tied, lexicographically smallest transcript_id.
        """
        tags = ['MANE_Select', 'CCDS',
                'appris_principal_1', 'Ensembl_canonical']
        df = df.copy()

        # 1) Compute total CDS length per transcript
        cds = df[df['feature'] == 'CDS'].copy()
        cds['cds_len'] = cds['end'] - cds['start'] + 1
        cds_len = (
            cds.groupby('transcript_id')['cds_len']
            .sum()
            .rename('total_cds_len')
        )

        # 2) Extract boolean flags for each priority tag
        df['tag'] = df['tag'].fillna('').astype(str)
        for t in tags:
            df[t] = df['tag'].str.contains(rf'\b{t}\b')

        # 3) Build one-row-per-transcript summary
        tx = (
            df[['gene_id', 'transcript_id', 'gene_name']]
            .drop_duplicates()
            .set_index('transcript_id')
        )
        tx = tx.join(cds_len, how='left').fillna({'total_cds_len': 0})
        flags = (
            df.drop_duplicates(subset=['transcript_id'])
            .set_index('transcript_id')[tags]
        )
        tx = tx.join(flags)

        # 4) Internal function to select the best transcript per gene
        def select_for_gene(group: pd.DataFrame) -> pd.Series:
            # a) Compute priority_score for each transcript
            def score_row(r):
                if r['MANE_Select']:
                    return 0
                if r['CCDS']:
                    return 1
                if r['appris_principal_1']:
                    return 2
                if r['Ensembl_canonical']:
                    return 3
                return 4

            group = group.copy()
            group['priority_score'] = group.apply(score_row, axis=1)

            # a) Filter to transcripts with the minimum priority_score
            min_score = group['priority_score'].min()
            candidates = group[group['priority_score'] == min_score].copy()

            # b) Hierarchical tag filtering
            for t in tags:
                has_t = candidates[t]
                # If this tag splits the set, keep only those that have it
                if has_t.any() and not has_t.all():
                    candidates = candidates[has_t]

            # c) If still tied, keep those with the highest CDS length
            if len(candidates) > 1:
                max_len = candidates['total_cds_len'].max()
                candidates = candidates[candidates['total_cds_len'] == max_len]

            # d) If still tied, pick the lexicographically smallest transcript_id
            selected = candidates.sort_index().iloc[0]

            # Add the transcript_id to the result
            selected = selected.copy()
            # The index is the transcript_id
            selected['transcript_id'] = selected.name
            return selected

        # Apply per-gene selection and return one row per gene
        canonical = (
            tx
            .reset_index()
            .groupby('gene_id')
            .apply(lambda g: select_for_gene(g.set_index('transcript_id')))
            .reset_index(drop=True)
        )

        ordered_cols = [
            'gene_id',
            'transcript_id',
            'gene_name',
            'total_cds_len',
            'MANE_Select',
            'CCDS',
            'appris_principal_1',
            'Ensembl_canonical',
            'priority_score'
        ]

        canonical = canonical[ordered_cols]

        return canonical


class SequenceExtractor:
    """Handles sequence extraction and processing, preserving original variable names and output schema."""

    def __init__(self, genome):
        self.genome = genome

    def fetch_sequence(self, chrom: str, start: int, end: int, strand: str, L: int) -> str:
        """
        Fetch chrom:start-end (1-based inclusive) with ±L padding.
        Returns the sequence string on the correct strand.
        """
        # pyfaidx is 1-based inclusive under the hood
        seq = self.genome[chrom][start - 1 - L: end + L].seq
        if strand == '-':
            comp = str.maketrans("ACGTacgt", "TGCAtgca")
            seq = seq.translate(comp)[::-1]
        return seq

    def add_exon_sequences(
        self,
        canonical_df: pd.DataFrame,
        annotation_df: pd.DataFrame,
        chromosome: str,
        L: int
    ) -> pd.DataFrame:
        """
        For each (gene_id, transcript_id) pair in `canonical_df`:
          1. Filter exons in `annotation_df` for the specified `chromosome`, `gene_id`, and `transcript_id`.
          2. Sort exons by `exon_number`.
          3. For each exon:
             a. Expand its genomic coordinates by ±padding bases to capture flanking context.
             b. Retrieve the padded DNA sequence via `fetch_sequence`.
             c. Strip off the padding to obtain the core exon sequence.
          4. Assemble and record, per transcript:
             - dna_coord_exon_seqs: list of [start−padding, end+padding] for each exon.
             - dna_exon_seqs: list of padded exon sequences (with flanking padding).
             - dna_abs_coord_exon: list of lists of genomic positions for each padded exon.
             - rna_exon_seqs: concatenated core exon sequences (padding removed).
             - rna_abs_coord_exon_all_seq: flat list of genomic positions for all core exons.
             - rna_coord_exon: list of lists of RNA‐relative positions per exon (1‐based).
             - rna_coord_all_seq: flat list of RNA‐relative positions for the entire transcript.


        Parameters
        ----------
        canonical_df : pd.DataFrame
            Must contain columns ['gene_id', 'transcript_id'] for canonical isoforms.
        annotation_df : pd.DataFrame
            Must contain columns ['seqname', 'gene_id', 'transcript_id',
            'exon_number', 'start', 'end', 'strand'] with exon annotations.
        chromosome : str
            Chromosome to filter on (e.g. "chr1").
        padding : int
            Number of bases to pad upstream and downstream of each exon.

        Returns
        -------
        pd.DataFrame
            A copy of `canonical_df` with the following nine added columns:
              - dna_coord_exon_seqs
              - dna_coord_exon_seqs_x0
              - dna_coord_exon_seqs_xf
              - dna_abs_coord_exon
              - dna_exon_seqs
              - rna_exon_seqs
              - rna_abs_coord_exon_all_seq
              - rna_coord_exon
              - rna_coord_all_seq
        """

        # --- containers to accumulate ---
        all_coords_dna = []                 # [[ [start-L, end+L], ... ], ...] per transcript
        # [[seq_exon1, seq_exon2, ...], ...] per transcript
        all_exon_seqs_dna = []
        all_concat_rna = []  # RNA exon sequence (core exons only)
        first_starts = []                   # dna_coord_exon_seqs_x0
        end_starts = []                     # dna_coord_exon_seqs_xf
        # [[list_of_positions_exon1, exon2, ...], ...] per transcript (padded)
        coord_exons_final = []

        coord_abs_exons_final_rna = []     # rna_abs_coord_exon_all_seq
        rna_coord_relative_exone = []       # rna_coord_exon
        rna_coord_relative_all_seq = []     # rna_coord_all_seq

        # Iterate over each canonical isoform
        for _, can_row in canonical_df.iterrows():
            gene_i = can_row['gene_id']
            transcript_i = can_row['transcript_id']

            # Filter annotation for this chromosome, gene, and transcript; sort by exon_number
            exons = (
                annotation_df[
                    (annotation_df['seqname'] == chromosome) &
                    (annotation_df['gene_id'] == gene_i) &
                    (annotation_df['transcript_id'] == transcript_i)
                ]
                .sort_values('exon_number',  ascending=True)
            )

            # Lists to collect per-exon data for this transcript

            # list of [start-L, end+L] per exon
            coords = []
            # list of padded exon sequences per exon
            seqs = []
            # list of lists of genomic positions per padded exon
            coord_exones = []
            # list of core exon sequences (padding removed) per exon
            seqs_rna_exons = []
            # list of lists of genomic positions per core exon (padding removed)
            coord_exones_rna = []

            # Loop through each exon record
            # Note: 'start' and 'end' are original exon coordinates (1-based inclusive)
            for _, exon in exons.iterrows():
                s, e = int(exon['start']), int(exon['end'])

                # Expand coordinates by ±L to include flanking regions for 5-mer context
                padded_start = s - L
                padded_end = e + L
                coords.append([padded_start, padded_end])

                # Fetch the DNA sequence for [s-L, e+L] on the specified strand
                seq_padded = self.fetch_sequence(
                    exon['seqname'], s, e, exon['strand'], L
                )
                seqs.append(seq_padded)

                # Build a list of all genomic positions for the padded region [padded_start .. padded_end]
                padded_positions = list(range(padded_start, padded_end + 1))
                coord_exones.append(padded_positions)

                # Extract the core exon sequence (remove L bases from both ends of seq_padded)
                # This yields the actual exon sequence without padding.
                core_seq = seq_padded[L:-L]
                seqs_rna_exons.append(core_seq)

                # Core exon positions (without padding) correspond to the central region
                core_positions = padded_positions[L:-L]
                coord_exones_rna.append(core_positions)

            # Flatten the list of core exon positions to get one list of all exon positions → one list
            merge_abs_coord_exon_rna = np.hstack(coord_exones_rna).tolist()

            # Build RNA-relative coordinates for each exon: 1-based positions within the transcript
            rna_style_lists = []
            counter = 1
            for positions in coord_exones_rna:
                length = len(positions)

                # Assign a consecutive block of RNA-relative positions for this exon
                rna_style_lists.append(list(range(counter, counter + length)))
                counter += length

            # Flatten all exon-level RNA-relative lists into a single list
            merge_rna_style_lists = np.hstack(rna_style_lists).tolist()

            # Record results for this transcript
            all_coords_dna.append(coords)
            all_exon_seqs_dna.append(seqs)

            # Record the first padded start and last padded end (for sorting later)
            first_starts.append(coords[0][0] if coords else pd.NA)
            end_starts  .append(coords[-1][-1] if coords else pd.NA)

            # Concatenate all core exon sequences to form the full transcript sequence
            all_concat_rna.append(''.join(seqs_rna_exons))

            # Store flat lists of core exon genomic positions
            coord_abs_exons_final_rna.append(merge_abs_coord_exon_rna)
            # Store padded exon genomic positions per exon
            coord_exons_final.append(coord_exones)

            # Store RNA-relative exon coordinates per exon and for the full transcript
            rna_coord_relative_exone.append(rna_style_lists)
            rna_coord_relative_all_seq.append(merge_rna_style_lists)

        # assemble output
        out = canonical_df.copy()
        out['dna_coord_exon_seqs'] = all_coords_dna
        out['dna_coord_exon_seqs_x0'] = pd.to_numeric(
            first_starts, errors='coerce')
        out['dna_coord_exon_seqs_xf'] = pd.to_numeric(
            end_starts, errors='coerce')
        out['dna_abs_coord_exon'] = coord_exons_final
        out['dna_exon_seqs'] = all_exon_seqs_dna
        out['rna_exon_seqs'] = all_concat_rna
        out['rna_abs_coord_exon_all_seq'] = coord_abs_exons_final_rna
        out['rna_coord_exon'] = rna_coord_relative_exone
        out['rna_coord_all_seq'] = rna_coord_relative_all_seq

        return out.sort_values('dna_coord_exon_seqs_x0', ascending=True).reset_index(drop=True)


class MutationAnalyzer:
    """Handles k-mer extraction and single-base mutations exactly as in your original loops."""

    NUCLEOTIDES = ['A', 'C', 'G', 'T']

    def __init__(self, k_mer_neighbors: int = 2, codon_table: dict = None):
        """
        Parameters
        ----------
        k_mer_neighbors : int
            Number of flanking bases on each side of the central base (so k-mer length = 2*k_mer_neighbors+1).
        codon_table : dict
            A dict mapping 3-letter codons (uppercase) to amino acids, e.g. {'ATG': 'M', ...}.
        """
        self.k_mer_neighbors = k_mer_neighbors
        self.codon_table = codon_table or {}

    def generate_kmer_df(self, sequences_df: pd.DataFrame) -> pd.DataFrame:
        """
        Reproduces your original `df_kmers = pd.DataFrame(records)` step.

        Input DF must have columns:
          - dna_exon_seqs            (list of padded exon strings)
          - dna_abs_coord_exon       (list of lists of genomic positions)
          - rna_coord_exon           (list of lists of RNA coords per exon)
          - rna_coord_all_seq        (flat list of all RNA coords)
          - rna_exon_seqs        (full concatenated RNA sequence)

        Returns a DataFrame with columns:
          ['gene_id','transcript_id','gene_name',
           'k_mer','central_base','codon','exon_index',
           'coord_abs','coord_rna','mod3']
        """
        records = []
        k = self.k_mer_neighbors

        for row in sequences_df.itertuples(index=False):
            # row.gene_id and row.transcript_id are the IDs for this row
            gene_id = row.gene_id
            transcript_id = row.transcript_id
            gene_name = getattr(row, 'gene_name', None)

            # row.exon_seqs is assumed to be a list of exon‐strings, e.g. ["ATGC…", "CGTA…", …]
            exons = row.dna_exon_seqs
            coords_exon_all = row.dna_abs_coord_exon

            coords_rna = row.rna_coord_exon
            coords_rna_all = row.rna_coord_all_seq
            all_sequence_rna = row.rna_exon_seqs

            # Loop over each exon in that list
            for j, exon in enumerate(exons):
                coords_exon = coords_exon_all[j]
                coords_exon_rna = coords_rna[j]

                # exon is just a Python string, so we can slice it as before
                for i in range(k, len(exon) - k):
                    k_mers = exon[i-k: i+k+1]

                    central = k_mers[k]
                    central_coord = coords_exon[i]

                    # RNA index within this exon
                    rna_index = i - k
                    central_coord_rna = coords_exon_rna[rna_index]
                    mod3 = central_coord_rna % 3

                    # slice out the codon from the full RNA sequence
                    idx = coords_rna_all.index(central_coord_rna)
                    if mod3 == 1:
                        codon = all_sequence_rna[idx: idx+3]
                    if mod3 == 2:
                        codon = all_sequence_rna[idx-1: idx+2]
                    if mod3 == 0:
                        codon = all_sequence_rna[idx-2: idx+1]
                    # print(k - k_mer_neighbors, k, k_mers, central, codon,
                    #       central_coord, central_coord_rna, mod3, gene_id, transcript_id, f"exon_index={j}")

                    records.append({
                        'gene_id': gene_id,
                        'transcript_id': transcript_id,
                        'gene_name': gene_name,
                        'k_mer': k_mers,
                        'central_base': central,
                        'codon': codon,
                        'exon_index': j,
                        'coord_abs': central_coord,
                        'coord_rna': central_coord_rna,
                        'mod3': mod3
                    })

        return pd.DataFrame(records)

    def generate_mutation_df(self, df_kmers: pd.DataFrame) -> pd.DataFrame:
        """
          - iterate df_kmers
          - build all single-base mutations
          - annotate wild vs mutant AA and is_nonsynonymous
        """
        new_records = []
        nts = self.NUCLEOTIDES

        for _, row in df_kmers.iterrows():
            wt_base = row['central_base']
            codon = row['codon']
            # Convert mod3 into a 0-based index within the codon:
            #   mod3 = 1 → index = 0
            #   mod3 = 2 → index = 1
            #   mod3 = 0 → index = 2
            frame_index = int(row['mod3']) - 1
            if frame_index < 0:
                frame_index = 2  # when mod3 == 0

            # Get the wild-type amino acid
            # wt_aa = self.codon_table.get(codon.upper(), None)
            wt_aa = CodonTable.translate(codon)

            for nt in nts:
                if nt == wt_base:
                    continue  # Skip the wild-type base

                # Build the mutated codon by replacing at the correct position
                codon_list = list(codon)
                codon_list[frame_index] = nt
                mutated_codon = ''.join(codon_list)

                # Get the mutant amino acid
                # mut_aa = self.codon_table.get(mutated_codon, None)
                mut_aa = CodonTable.translate(mutated_codon)

                # Determine if the mutation is synonymous (0) or nonsynonymous (1)
                # If either codon is not in the table, treat as nonsynonymous
                is_nonsyn = 0 if (
                    wt_aa is not None and mut_aa is not None and wt_aa == mut_aa) else 1

                # Determine if the mutation is a nonsynonymous stop codon
                is_nonsyn_stop = int(is_nonsyn and mut_aa != '*')

                # Copy all original fields and add the mutation
                rec = row.to_dict()
                rec.update({
                    'mutated_base': nt,
                    'mutated_codon': mutated_codon,
                    'wild_aa': wt_aa,
                    'mutated_aa': mut_aa,
                    'is_nonsynonymous': is_nonsyn,
                    'is_nonsynonymous_nonstop': is_nonsyn_stop
                })
                new_records.append(rec)

        # enforce desired column order:
        desired_cols = [
            'gene_id', 'transcript_id', 'gene_name',
            'k_mer', 'central_base', 'mutated_base',
            'codon', 'mutated_codon', 'wild_aa', 'mutated_aa', 'is_nonsynonymous', 'is_nonsynonymous_nonstop',
            'exon_index', 'coord_abs', 'coord_rna', 'mod3'
        ]
        df_mut = pd.DataFrame(new_records)
        return df_mut[desired_cols]


class PeptideAnalyzer:
    """
    Generates mutated peptides from RNA sequences and a mutations DataFrame,
    following the structure of your original script.
    """

    def __init__(self, peptide_length: int = 9):
        """
        Parameters
        ----------
        peptide_length : int
            Number of amino acids per peptide.
        """
        self.peptide_length = peptide_length
        self.window_length = 3 * peptide_length

    def generate_sequence_vectors(
        self,
        sequences: List[str]
    ) -> Tuple[List[List[int]], List[List[int]]]:
        """
        For each RNA sequence, build
          1. a list of consecutive 1-based positions (vec_num),
          2. a list of codon-grouped positions (vec_grouped: 1,1,1,2,2,2,…).

        Returns
        -------
        vec_num : List[List[int]]
        vec_grouped : List[List[int]]
        """
        vec_num = []
        counter = 1
        for seq in sequences:
            length = len(seq)
            vec_num.append(list(range(counter, counter + length)))
            counter += length

        # build a flat codon‐grouped vector
        total_length = sum(len(seq) for seq in sequences)
        num_blocks = (total_length + 2) // 3
        flat_group = []
        for block in range(1, num_blocks + 1):
            flat_group.extend([block] * 3)

        vec_grouped = []
        start = 0
        for seq in sequences:
            end = start + len(seq)
            vec_grouped.append(flat_group[start:end])
            start = end

        return vec_num, vec_grouped

    def generate_mutated_peptides(
        self,
        sequences_df: pd.DataFrame,
        mutations_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        1) Flatten sequences_df into lists of:
             - gene_ids_t
             - transcript_ids_t
             - rna_seqs_all_t
             - abs_coord_all_t
        2) Build RNA coordinate vectors.
        3) Slide a window of self.window_length over each RNA, stepping by 3.
        4) For each peptide window, scan through df_mutated for matching coords,
           generate the mutated peptide, and collect metadata.

        Input DF must contain:
          - gene_id, transcript_id
          - rna_exon_seqs               (full CDS RNA)
          - rna_abs_coord_exon_all_seq      (flat genomic coords)
        mutations_df must contain:
          - gene_id, transcript_id, coord_abs, mutated_base, central_base,
            codon, mutated_codon, wild_aa, mutated_aa, is_nonsynonymous

        Returns
        -------
        pd.DataFrame
            One row per mutated peptide, with columns:
            ['gene_id','transcript_id','gene_name',
             'original_peptide','mutated_peptide',
             'original_peptide_aa','mutated_peptide_aa',
             'codon_position','base_position_in_codon',
             'absolute_coord','wild_codon','mutated_codon',
             'wild_aa','mutated_aa','is_nonsynonymous',
             'position_in_peptide',
             'coord_seq_peptide_x0','coord_seq_peptide_xf']
        """
        # 1) collect per-transcript lists
        gene_ids_t: List[List[str]] = []
        transcript_ids_t: List[List[str]] = []
        rna_seqs_all_t: List[str] = []
        abs_coord_all_t: List[List[int]] = []

        for row in sequences_df.itertuples(index=False):
            seq = row.rna_exon_seqs
            coords = row.rna_abs_coord_exon_all_seq
            length = len(seq)

            gene_ids_t.append([row.gene_id] * length)
            transcript_ids_t.append([row.transcript_id] * length)
            rna_seqs_all_t.append(seq.upper())
            abs_coord_all_t.append(coords)

        # 2) build RNA coordinate vectors
        rna_n_coordinates, _ = self.generate_sequence_vectors(
            rna_seqs_all_t)

        peptides_mut: List[Dict] = []

        # 3) slide window over each sequence

        # for idx, seq in enumerate(rna_seqs_all_t):
        # for idx, seq in tqdm(enumerate(rna_seqs_all_t), total=len(rna_seqs_all_t), desc="processing transcripts"):

        for idx, seq in tqdm(enumerate(rna_seqs_all_t), total=len(rna_seqs_all_t),
                             desc="processing transcripts", file=sys.stdout):
            seq = seq.upper()
            # seq = dna_seq.upper()[0:42]  # Take the first 30 bases as an example (same as before)
            # Take the first 30 bases as an example (same as before)
            # seq = dna_seq.upper()[0:30]

            rna_coords = rna_n_coordinates[idx]
            abs_coords = abs_coord_all_t[idx]
            gene_id = gene_ids_t[idx][0]
            transcript_id = transcript_ids_t[idx][0]

            # filter by gene, transcript
            filter_mutations = (
                (mutations_df['gene_id'] == gene_id) &
                (mutations_df['transcript_id'] == transcript_id))
            gene_mutation = mutations_df[filter_mutations]

            # slide window
            for start in range(0, len(seq) - self.window_length + 1, 3):
                peptide = seq[start:start + self.window_length]
                pep_coords = abs_coords[start:start + self.window_length]
                pep_coords_x0 = pep_coords[0]
                pep_coords_xf = pep_coords[-1]

                # 4) scan for mutations within this window
                # Filtering mutation by absolute coordinate
                """
                window_mutations = mutations_df[
                    (mutations_df['gene_id'] == gene_id) &
                    (mutations_df['transcript_id'] == transcript_id) &
                    (mutations_df['coord_abs'].isin(pep_coords))
                 ]
                """
                window_mutations = gene_mutation[
                    (gene_mutation['coord_abs'].isin(pep_coords))
                ]

                for mut in window_mutations.itertuples(index=False):
                    mut_pos = pep_coords.index(mut.coord_abs)
                    codon_num = mut_pos // 3
                    pos_in_codon = mut_pos % 3

                    # build mutated peptide
                    mutated_peptide = (
                        peptide[:mut_pos] +
                        mut.mutated_base +
                        peptide[mut_pos + 1:]
                    )

                    peptides_mut.append({
                        'gene_id': gene_id,
                        'transcript_id': transcript_id,
                        'gene_name': mut.gene_name,
                        'wild_peptide': peptide,
                        'mutated_peptide': mutated_peptide,
                        'wild_peptide_aa': CodonTable.translate_sequence(peptide),
                        'mutated_peptide_aa': CodonTable.translate_sequence(mutated_peptide),
                        'codon_position': codon_num,
                        'base_position_in_codon': pos_in_codon,
                        'absolute_coord': mut.coord_abs,
                        'wild_codon': mut.codon,
                        'mutated_codon': mut.mutated_codon,
                        'wild_aa': mut.wild_aa,
                        'mutated_aa': mut.mutated_aa,
                        'is_nonsynonymous': mut.is_nonsynonymous,
                        'is_nonsynonymous_nonstop': mut.is_nonsynonymous_nonstop,
                        'position_in_peptide': mut_pos,
                        'coord_seq_peptide_x0': pep_coords_x0,
                        'coord_seq_peptide_xf': pep_coords_xf
                        # 'rna_n_coordinates': rna_coords
                    })

        df_all_peptides = pd.DataFrame(peptides_mut)

        # Select unique mutated peptides whose amino acid sequences are nonsynonymous and not stop codons
        df_unique_nonsyn_nonstop = (
            df_all_peptides[df_all_peptides["is_nonsynonymous_nonstop"] == 1]
            .drop_duplicates(subset=["mutated_peptide_aa"])
            .loc[:, ["wild_peptide_aa", "mutated_peptide_aa"]]
            .reset_index(drop=True)
        )

        return df_all_peptides, df_unique_nonsyn_nonstop


class GenomicsAnalysisPipeline:
    """Main pipeline class that orchestrates the analysis."""

    def __init__(self, config: GenomicsConfig):
        self.config = config
        self.genome = None
        self.annotations = None

        # Initialize components
        self.isoform_selector = IsoformSelector()
        self.mutation_analyzer = MutationAnalyzer(config.k_mer_neighbors)
        self.peptide_analyzer = PeptideAnalyzer(config.peptide_length)

    def load_data(self):
        """Load genome and annotation data."""
        print("Loading genome and annotation data...")

        # Load genome
        self.genome = Fasta(self.config.path_grch38)
        self.sequence_extractor = SequenceExtractor(self.genome)

        # Load annotations
        self.annotations = read_gtf(
            self.config.path_cds_annotation).to_pandas()
        self.annotations["transcript_support_level"] = (
            self.annotations["transcript_support_level"]
            .replace(r'^\s*$', pd.NA, regex=True)
        )

        print(f"Loaded {len(self.annotations)} annotation records")

    def run_analysis(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Run the complete analysis pipeline."""
        if self.genome is None or self.annotations is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        print(f"Running analysis for chromosome {self.config.chromosome}")

        # Step 1: Filter annotations and select canonical isoforms
        print("Step 1: Selecting canonical isoforms...")
        chr_i = self.config.chromosome

        filtered_annotations = self.isoform_selector.filter_isoforms_by_phase(
            self.annotations, self.config.chromosome
        )

        canonical_transcripts = self.isoform_selector.pick_canonical_transcript(
            filtered_annotations
        )
        print(f"Selected {len(canonical_transcripts)} canonical transcripts")

        # print('annotations', self.annotations)
        # print('canonical_transcripts', canonical_transcripts)

        # Step 1-b: compute your two gene→isoform summaries
        summary_all_annot = IsoformSelector.get_isoforms_by_gene(
            self.annotations, chr_i
        )
        summary_filter_annot = IsoformSelector.get_isoforms_by_gene(
            filtered_annotations, chr_i
        )
        # print('summary_all_annot', summary_all_annot)
        # print('summary_filter_annot', summary_filter_annot)

        # Step 2: Extract sequences
        print("Step 2: Extracting sequences...")
        sequences_with_data = self.sequence_extractor.add_exon_sequences(
            canonical_transcripts,
            filtered_annotations,
            self.config.chromosome,
            self.config.k_mer_neighbors
        ).head(2)
        print(
            f"Extracted sequences for {len(sequences_with_data)} transcripts")

        # Step 3: Generate k-mers and mutations
        print("Step 3a: Generating k-mer table…")
        df_kmers = self.mutation_analyzer.generate_kmer_df(
            sequences_with_data)
        print(f"→ {len(df_kmers)} k-mers generated")

        print("Step 3b: Generating mutation table…")
        mutations_df = self.mutation_analyzer.generate_mutation_df(
            df_kmers
        )
        print(f"Generated {len(mutations_df)} mutations")

        # Step 4: Generate peptides
        print("Step 4: Generating mutated peptides...")
        peptides_df, df_unique_nonsyn_nonstop = self.peptide_analyzer.generate_mutated_peptides(
            sequences_with_data, mutations_df
        )
        print(f"Generated {len(peptides_df)} mutated peptides")

        return sequences_with_data, df_kmers, mutations_df, peptides_df, df_unique_nonsyn_nonstop, summary_all_annot, summary_filter_annot

    def get_summary_statistics(self, sequences_df: pd.DataFrame,
                               mutations_df: pd.DataFrame,
                               peptides_df: pd.DataFrame) -> Dict:
        """Generate summary statistics."""
        stats = {
            'n_genes': sequences_df['gene_id'].nunique(),
            'n_transcripts': sequences_df['transcript_id'].nunique(),
            'n_mutations': len(mutations_df),
            'n_nonsynonymous_mutations': len(mutations_df[mutations_df['is_nonsynonymous'] == 1]),
            'n_peptides': len(peptides_df),
            'n_nonsynonymous_peptides': len(peptides_df[peptides_df['is_nonsynonymous'] == 1])
        }
        return stats


# %%
"""
config = GenomicsConfig()
genome = Fasta(config.path_grch38)
se = SequenceExtractor(genome)
chrom = config.chromosome

df_annotations = read_gtf(config.path_cds_annotation).to_pandas()
# (opcional clean-up transcript_support_level si lo usas)
if "transcript_support_level" in df_annotations:
    df_annotations["transcript_support_level"] = (
        df_annotations["transcript_support_level"]
        .replace(r'^\s*$', pd.NA, regex=True)
    )

# 2) Filtrar isoformas por fase y escoger transcript canónico (igual que antes)

df_valid_annotations = IsoformSelector.filter_isoforms_by_phase(
    df_annotations, config.chromosome)

df_valid_annotations[(df_valid_annotations['seqname'] == "chr1") &
                     (df_valid_annotations['gene_id'] == "ENSG00000000457.14") &
                     (df_valid_annotations['transcript_id'] == "ENST00000367771.11")]



canonical = IsoformSelector.pick_canonical_transcript(df_valid_annotations)


canonical_with_exons_penta = se.add_exon_sequences(
    canonical_df=canonical,
    annotation_df=df_valid_annotations,
    chromosome=config.chromosome,
    L=config.k_mer_neighbors
)





canonical_with_exons_penta = se.add_exon_sequences(
    canonical_df=canonical,
    annotation_df=df_valid_annotations,
    chromosome=config.chromosome,
    L=config.k_mer_neighbors
)


# 3) Initialize the mutation analyzer (only needs k_mer_neighbors)
ma = MutationAnalyzer(
    k_mer_neighbors=config.k_mer_neighbors
)

# 4) First, build the k-mer DataFrame exactly as before:
df_kmers_1 = ma.generate_kmer_df(canonical_with_exons_penta.head(4))

# 5) Then derive all single-base mutations from it:
df_mutated_1 = ma.generate_mutation_df(df_kmers_1)

pa = PeptideAnalyzer(peptide_length=config.peptide_length)

# Esto mostrará un tqdm por cada transcript procesado:
df_peptides_1 = pa.generate_mutated_peptides(
    sequences_df=canonical_with_exons_penta.head(4),
    mutations_df=df_mutated_1
)

df_peptides_33 = df_peptides[
    ['gene_id', 'transcript_id', 'gene_name', 'wild_peptide',
     'mutated_peptide', 'wild_peptide_aa', 'mutated_peptide_aa',
     'codon_position', 'base_position_in_codon', 'absolute_coord',
     'wild_codon', 'mutated_codon', 'wild_aa', 'mutated_aa',
     'is_nonsynonymous', 'position_in_peptide', 'coord_seq_peptide_x0',
     'coord_seq_peptide_xf']
]
"""

# %%
# 1) Construye la configuración con valores por defecto
config = GenomicsConfig()
print("=== CONFIGURATION ===")
print(f" GTF file  : {config.path_cds_annotation}")
print(f" FASTA file: {config.path_grch38}")
print(f" Chromosome: {config.chromosome}")
print(f" k-mer pad : {config.k_mer_neighbors}")
print(f" Peptide length: {config.peptide_length}")
print("-" * 40)

# 2) Instancia y carga datos
pipe = GenomicsAnalysisPipeline(config)
pipe.load_data()

# 3) Ejecuta el pipeline completo
(canonical_with_exons_penta,
    df_kmers,
    df_mutated,
    df_peptides,
    df_unique_nonsyn_nonstop,
    summary_all_annot,
    summary_filter_annot) = pipe.run_analysis()

# 4) Calcula estadísticas y muéstralas
stats = pipe.get_summary_statistics(
    canonical_with_exons_penta,
    df_mutated,
    df_peptides
)
print("\n=== SUMMARY ===")
for name, count in stats.items():
    print(f"{name:30s}: {count}")

# 5) Opcional: mira los primeros registros de cada DataFrame
print("\n--- sample outputs ---")
print("sequences (first 2 rows):")
print(canonical_with_exons_penta, "\n")

print("kmers (first 5 rows):")
print(df_kmers.head(5), "\n")

print("mutations (first 5 rows):")
print(df_mutated.head(5), "\n")

print("peptides (first 5 rows):")
print(df_peptides.head(5), "\n")

# %%
# %%
# %%
# %%
# %%
# save outputs
# save outputs in Parquet
outdir = Path(
    '/Users/fordonez/Documents/PhD_Thesis/somatic_mutation_generation/results_mutations')
outdir.mkdir(parents=True, exist_ok=True)

canonical_with_exons_penta.to_parquet(
    outdir/f"{'chr1'}_sequences.parquet", index=False
)
df_kmers.to_parquet(
    outdir/f"{'chr1'}_kmers.parquet", index=False
)
df_mutated.to_parquet(
    outdir/f"{'chr1'}_mutations.parquet", index=False
)
df_peptides.to_parquet(
    outdir/f"{'chr1'}_peptides.parquet", index=False
)
# %%
# %%
# a = 10
# df_peptides[["is_nonsynonymous"]][a+0:80+a].sum()


'''

def main():
    parser = argparse.ArgumentParser(
        description="Genomics pipeline: isoforms → k-mers → peptides"
    )
    parser.add_argument("--gtf",    required=True, help="CDS GTF file")
    parser.add_argument("--fasta",  required=True, help="Genome FASTA file")
    parser.add_argument("--chrom",  default="chr1",
                        help="Chromosome to analyse")
    parser.add_argument("--kmer",   type=int, default=2,
                        help="k-mer neighbor padding")
    parser.add_argument("--pep-len", type=int, default=9,
                        help="Peptide length (AAs)")
    parser.add_argument("--outdir",   default="results",
                        help="Directory in which to write output Parquet files")
    args = parser.parse_args()

    # build & override config
    config = GenomicsConfig()
    config.path_cds_annotation = args.gtf
    config.path_grch38 = args.fasta
    config.chromosome = args.chrom
    config.k_mer_neighbors = args.kmer
    config.peptide_length = args.pep_len
    config.window_length = 3 * config.peptide_length

    # run
    pipe = GenomicsAnalysisPipeline(config)
    pipe.load_data()

    (canonical_with_exons_penta,
     df_kmers,
     df_mutated,
     df_peptides, 
     summary_all_annot,
     summary_filter_annot) = pipe.run_analysis()

    stats = pipe.get_summary_statistics(
        canonical_with_exons_penta, df_mutated, df_peptides
    )
    print("\n=== SUMMARY ===")
    for k, v in stats.items():
        print(f"{k:30s}: {v}")

    # save outputs
    # save outputs in Parquet
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    canonical_with_exons_penta.to_parquet(
        outdir/f"{args.chrom}_sequences_table.parquet", index=False
    )
    df_kmers.to_parquet(
        outdir/f"{args.chrom}_kmers.parquet", index=False
    )
    df_mutated.to_parquet(
        outdir/f"{args.chrom}_mutations.parquet", index=False
    )
    df_peptides.to_parquet(
        outdir/f"{args.chrom}_peptides.parquet", index=False
    )
    """
    outdir = Path("results")
    outdir.mkdir(exist_ok=True)
    canonical_with_exons_penta.to_csv(outdir/"sequences.csv", index=False)
    df_kmers   .to_csv(outdir/f"{args.chrom}_kmers.csv",      index=False)
    df_mutated .to_csv(outdir/f"{args.chrom}_mutations.csv",  index=False)
    df_peptides.to_csv(outdir/f"{args.chrom}_peptides.csv",   index=False)
    """
    print(f"\nAll results written to {outdir.resolve()}")

if __name__ == "__main__":
    main()


python pipeline.py \
  --gtf /path/to/gencode.v47.CDS.gtf \
  --fasta /path/to/GRCh38.p14.genome.fa \
  --chrom chr1 \
  --kmer 2 \
  --pep-len 9

'''

# %%
