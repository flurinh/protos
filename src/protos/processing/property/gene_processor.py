from protos.loaders.uniprot_loader import *
import os
import pandas as pd


class GeneProcessor:
    def __init__(self,
                 gene_table='gene_table_V1.csv',
                 mapping_table='mapping_table.csv',
                 path='data/gene/',
                 preload=True):
        self.path = path
        self.gene_table = gene_table
        self.mapping_table = mapping_table
        self.data = pd.DataFrame()
        self.mapping = pd.DataFrame()
        if preload:
            self.load_gene_table()

    def load_gene_table(self, load_mapping=True):
        filename_gene_table = os.path.join(self.path, self.gene_table)
        self.data = pd.read_csv(filename_gene_table)
        if load_mapping:
            filename_mapping_table = os.path.join(self.path, self.mapping_table)
            self.mapping = pd.read_csv(filename_mapping_table)

    def __len__(self):
        return len(self.data)

    def get_entries_by_uid(self, query_id):
        mask = self.data['uid'].str.contains(query_id, case=False, na=False)
        return self.data[mask]

    def get_uid_pdb_mapping(self):
        job_id = submit_id_mapping(
            from_db="UniProtKB_AC-ID", to_db="PDB", ids=self.data.uid.tolist()
        )
        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link + "?format=tsv")
            self.mapping = pd.DataFrame([(x.split('\t')[0], (x.split('\t')[1])) for x in results][1:], columns=['uid', 'pdb_id'])
