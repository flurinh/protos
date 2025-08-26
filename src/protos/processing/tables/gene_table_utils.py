import pandas as pd
from tqdm import trange
from protos.processing.grn.grn_table_utils import *
import time

# This should only be called once...

def create_default_gene_table(gene_table_cols, fam_df, limit = 10):
    gene_table_data = []
    for row in trange(len(fam_df)):
        if row < limit:
            gene_data = fam_df.iloc[row].to_dict()
            # get uniprot data (in this case just by getting the uniprot ID)
            uniprot_id = gene_data['uniprot']
            uniprot_query = []
            if uniprot_id == '':
                gene = gene_data['gene']
                try:
                    uniprot_query = get_uniprot(gene + gene_data['full_name'], batchsize=1)
                    if len(uniprot_query) > 0:
                        uniprot_query = uniprot_query.iloc[0].to_dict()
                    else:
                        print("error", gene + gene_data['full_name'])
                except:
                    print("error with link")
            else:
                try:
                    gene = gene_data['gene']
                    uniprot_query = get_uniprot(uniprot_id, batchsize=1)
                    if len(uniprot_query) > 0:
                        uniprot_query = uniprot_query.iloc[0].to_dict()
                        uniprot_id = uniprot_query['Entry']
                    else:
                        print("error", uniprot_id)
                except:
                    print("error with link")
            if len(uniprot_query) > 0:
                entity = uniprot_query['Entry Name']
                length = uniprot_query['Length']
                seq = uniprot_query['Sequence']
                organism = uniprot_query['Organism']
                f1 = gene_data['f1']
                f1_id = gene_data['family_id_0']
                f2 = gene_data['f2']
                f2_id = gene_data['family_id_1']
                f3 = gene_data['f3']
                f3_id = gene_data['family_id_2']
                other_names = gene_data['full_name']
                gene_table_row = [entity, gene, organism, length, seq,
                                  f1, f2, f3, f1_id, f2_id, f3_id,
                                  uniprot_id, '', '', '', other_names, []]
                gene_table_data.append(gene_table_row)
            time.sleep(1)
    return pd.DataFrame(gene_table_data, columns=gene_table_cols)

