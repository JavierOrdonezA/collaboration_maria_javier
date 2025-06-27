import numpy as np
import pandas as pd


class NeutralEvolutionSimulator:
    def __init__(self, n_subs, coord_file, cbase_output, syn_probs, non_syn_probs, output):
        self.Nsubs = n_subs
        self.coord_file = coord_file
        self.cbase_output = cbase_output
        self.syn_probs = syn_probs
        self.non_syn_probs = non_syn_probs
        self.output = output

        # Constants
        self.bases = ['A', 'C', 'G', 'T']
        self.trinucleotides = [a + b + c for a in self.bases for b in self.bases for c in self.bases]

        # Load and preprocess data
        self.load_data()

    def load_data(self):
        self.gene_annotations_df = pd.read_csv(self.coord_file, index_col=0)
        self.cbase_output_df = pd.read_csv(self.cbase_output, index_col=0)

        expressed_genes = self.gene_annotations_df['gene'].unique()
        print(len(expressed_genes))
        self.cbase_output_df = self.cbase_output_df[self.cbase_output_df.index.isin(expressed_genes)]
        self.list_of_genes = self.cbase_output_df.index.unique().tolist()
        print(len(self.list_of_genes))
        self.syn_df = pd.read_csv(self.syn_probs, names=['gene', 'pos', 'alt', 'prob'], dtype={"prob": "float64"})
        self.non_syn_df = pd.read_csv(self.non_syn_probs, names=['gene', 'pos', 'alt', 'prob'], dtype={"prob": "float64"})

        self.syn_by_gene = dict(tuple(self.syn_df.groupby('gene')))
        self.non_syn_by_gene = dict(tuple(self.non_syn_df.groupby('gene')))
        self.annotations_by_gene = dict(tuple(self.gene_annotations_df.groupby('gene')))

    def find_across_gene_probs(self):
        probs = self.cbase_output_df[['lambda_s', 'dN_lambda_n']].to_numpy().flatten()
        total = probs.sum()
        return (probs / total).tolist()

    #need to check if this keeps correct order
    def distribute_muts_across_genes(self, genes_probs):
        mutations = np.random.multinomial(self.Nsubs, genes_probs)
        #from index 0 in steps of 2
        muts_synon = mutations[::2]
        #from index 1 in steps of 2
        muts_non_synon = mutations[1::2]
        data = {'gene': self.list_of_genes, 'synon_muts': muts_synon, 'non_synon_muts': muts_non_synon}
        return pd.DataFrame(data).set_index('gene')

    def apply_muts_per_gene(self, gene, sim_df, gene_df):
        results = []
        alt_map = dict(zip(sim_df['pos'], sim_df['alt']))
        for _, exon in gene_df.iterrows():
            exon_start = exon['start']
            exon_seq = list(exon['seq'])
            mutated_seq = ''.join(alt_map.get(exon_start + i, base) for i, base in enumerate(exon_seq))
            results.append({
                'chr': exon['chr'],
                'start': exon_start,
                'end': exon['end'],
                'gene': gene,
                'strand': exon['strand'],
                'simulated_seq': mutated_seq
            })
        results_df = pd.DataFrame(results)
        results_df.to_csv(self.output, mode='a', index=False, header=False)

    def run(self):
        genes_probs = self.find_across_gene_probs()
        muts_per_gene_df = self.distribute_muts_across_genes(genes_probs)
        all_mutated_genes = set()

        for gene, row in muts_per_gene_df.iterrows():
            syn_gene_df = self.syn_by_gene.get(gene)
            non_syn_gene_df = self.non_syn_by_gene.get(gene)
            gene_df = self.annotations_by_gene.get(gene)
            if syn_gene_df is None or non_syn_gene_df is None or gene_df is None:
                continue

            Ngene_syn = row['synon_muts']
            Ngene_non_syn = row['non_synon_muts']

            syn_gene_df = syn_gene_df.copy()
            non_syn_gene_df = non_syn_gene_df.copy()

            syn_gene_df['sim_mut'] = np.random.multinomial(Ngene_syn, syn_gene_df['prob'].values)
            non_syn_gene_df['sim_mut'] = np.random.multinomial(Ngene_non_syn, non_syn_gene_df['prob'].values)

            gene_simulated_subs_df = pd.concat([syn_gene_df, non_syn_gene_df]).sort_index()
            gene_no0_simulated_subs_df = gene_simulated_subs_df[gene_simulated_subs_df.sim_mut != 0]

            if not gene_no0_simulated_subs_df.empty:
                all_mutated_genes.add(gene)
                self.apply_muts_per_gene(gene, gene_no0_simulated_subs_df, gene_df)

        # Save non-mutated genes
        genes_no_mut = [g for g in self.list_of_genes if g not in all_mutated_genes]
        no_mut_df = self.gene_annotations_df[self.gene_annotations_df['gene'].isin(genes_no_mut)]
        no_mut_df.to_csv(self.output, mode='a', index=False, header=False)


