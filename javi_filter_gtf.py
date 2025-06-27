import pandas as pd
from gtfparse import read_gtf


path_cds_annotation ='/home/maria/data/gencode.v47.annotation.gtf.gz'
extract_features = ['CDS','stop_codon']
output = '/home/maria/filter_transcripts/output/exon.bed'

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
        print(len(is_pc))

        # pass_pc = is_pc & df_chr['CCDS']
        pass_pc = is_pc & (df_chr['CCDS'] | df_chr['MANE_Select'])
        print(len(pass_pc))

        #i do not want to keep these non protein coding genes
        # pass_npc = ~is_pc & df_chr['level'].isin([1, 2])
        pass_npc = is_pc & df_chr['level'].isin([1, 2])
        df_q = df_chr[pass_pc | pass_npc].copy()
        df_q = df_chr[pass_pc].copy()

        #MARIA ADD
        #think need to filter by feature, so summing the correct items: 
        df_q = df_q[df_q['feature'].isin(['CDS','stop_codon'])]

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

        
        print(len(df_valid))
        # drop the helper column        df_valid = df_valid.drop(columns=['segment_length'])

        df_valid["exon_number"] = pd.to_numeric(
            df_valid["exon_number"], errors="coerce")
        #i am also oututting df_q to have the option to extract coordinates of other features ie introns, 
        #df_valid is now filtered only for ccds plus stop codon
        return df_valid, df_q

    @staticmethod

    #not used anymore
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
    
class GenomicsAnalysisPipeline:
    """Main pipeline class that orchestrates the analysis."""

    def __init__(self, path_cds_annotation, extract_features, output):
        self.path_cds_annotation = path_cds_annotation
        self.feature = extract_features
        self.output = output
        self.annotations = None

        # Initialize components
        self.isoform_selector = IsoformSelector()
        

    def load_data(self):
        """Load genome and annotation data."""
        print("Loading genome and annotation data...")

        
        # Load annotations
        #self.annotations = read_gtf(
            #self.path_cds_annotation).to_pandas()
        self.annotations = read_gtf(
            self.path_cds_annotation)
        self.annotations["transcript_support_level"] = (
            self.annotations["transcript_support_level"]
            .replace(r'^\s*$', pd.NA, regex=True)
        )

        print(f"Loaded {len(self.annotations)} annotation records")

    #my function to create bedfile of coordinates
    def create_bed_file(self, gtf_df):
        mask = gtf_df.apply(
        lambda row: self.transcript_gene_ids.get(row['gene_id']) == row['transcript_id'],
        axis=1
        )

        # Apply the mask
        transcript_df = gtf_df[mask]

        #filer based on feature
        feature_df =transcript_df[transcript_df['feature'].isin(self.feature)]

        #transcript_df = df[df['gene_id']==key and df['transcript_id']==value for key, value in transcript_gene_ids.]
        #filter also based on the transcript


        cut_df_exon = feature_df[['seqname','start', 'end', 'gene_name', 'strand']]
        #convert from gtf 1 based to bedfile 0
        cut_df_exon['start'] -= 1
        cut_df_exon.set_index('seqname', inplace=True)
        print(len(cut_df_exon['gene_name'].unique()))
        cut_df_exon.to_csv(self.output, header=False, sep ='\t')


    def run_analysis(self):
        """Run the complete analysis pipeline."""
        self.load_data()
        if self.annotations is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        

        # Step 1: Filter annotations and select canonical isoforms
        print("Step 1: Selecting canonical isoforms...")
        chr_i = 'chr1'

        #want to check how many isoforms i am losing here
        filtered_annotations, filtered_all_features = self.isoform_selector.filter_isoforms_by_phase(
            self.annotations, chr_i
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

        #creating dict of selected transcript id per gene
        self.transcript_gene_ids = canonical_transcripts.set_index('gene_id')['transcript_id'].to_dict()
        
        #intput df with protein coding filter, to select transcript and feature from, then save
        self.create_bed_file(filtered_all_features)
        return canonical_transcripts[['gene_id', 'transcript_id']]

transcripts=GenomicsAnalysisPipeline(path_cds_annotation, extract_features, output)
canoncial_transcripts = transcripts.run_analysis()
print(canoncial_transcripts)