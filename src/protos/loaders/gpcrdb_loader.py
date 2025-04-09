from data_download.utils import *
from data_download.gpcrdb.gpcrdb_loader_utils import *

from tqdm import tqdm
import pandas as pd
from pathlib import Path


COLS = ['cross', 'filler', 'uniprot(gene)', 'filler2', 'receptor family', 'Cl.', 'Species', 'Method', 'PDB', \
        'Refined Structure', 'Resolution', 'Preferred Chain', 'State', 'Degree active %', '% of Seq1', 'Family', \
        'Subtype', 'Note', '% of Seq2', 'Fusion', 'Antibodies', 'Name1', 'Type1', 'Function', 'Name2', 'Type2', \
        'D2x50 - S3x39', 'Sodium in structure', 'Authors', 'Reference', 'PDB date', 'Annotated']

regions = ['TM1', 'ICL1', 'TM2', 'ECL1', 'TM3', 'ICL2', 'TM4', 'ECL2',  'TM5', 'TM6', 'TM7', 'H8']
# gproteins: https://github.com/protwis/gpcrdb_data/blob/master/structure_data/annotation/g_proteins.csv


class GPCRdbDL:
    def __init__(self,
                 path='data/',
                 fileformat='pdb'):
        if fileformat == 'pdb':
            self.path_pdb = path + 'pdb/'
        elif fileformat == 'cif':
            self.path_pdb = path + 'mmcif/'
        else:
            raise ValueError("unknown format")

        if not Path(self.path_pdb).exists():
            raise OSError(f"Directory {self.path_pdb} not found")

        self.path_alignment = path + 'gpcrdb/'
        if not Path(self.path_alignment).exists():
            raise OSError(f"Directory {self.path_alignment} not found")

        self.path_table = path + 'gpcrdb/'

        self.fileformat = fileformat

        self.table = None
        self.filenames, self.pdb_ids = get_pdb_files(path=self.path_pdb)

    # ======================================================================================================================
    
    def download_pdbs(self, reload=False, update=False):
        pdb_table_ids = self.table['PDB'].tolist()
        missing = [x for x in pdb_table_ids if x not in self.pdb_ids]
        if reload or (len(missing) > 0):
            print("Reloading pdb files...")
            print("Missing mmcif_ret_microbial_opsins:", missing)
            for pdb in tqdm(missing):
                url = get_rcsb_download(pdb, self.fileformat)
                download_pdb(url, folder=self.path_pdb, fileformat=self.fileformat)
        elif update:
            self.update_pdbs()
        self.filenames, self.pdb_ids = get_pdb_files(path=self.path_pdb)

    def download_table(self, reload=True, filename='data_table.pkl'):
        table_path = self.path_table + filename
        if reload or (not os.path.isfile(self.path_table + filename)):
            self.table = get_table(reload=True, uniprot=False)
            self.table = self.table.drop(columns=['filler', 'filler2',
                                                  'Refined Structure',
                                                  '% of Seq1', 'Name2', 'Fusion', 'Note',
                                                  '% of Seq2', 'Antibodies', 'Name1', 'Type1', 'Type2',
                                                  'D2x50 - S3x39', 'Sodium in structure', 'Authors', 'Reference',
                                                  'PDB date', 'Annotated', 'pdb_link'])
        else:
            self.table = pd.read_pickle(self.path_table)
        print("writing gpcrdb table to file:", table_path)
        self.table.to_pickle(table_path)

    def add_row_table(self,
                      uniprot='',
                      receptor_family='',
                      cl='',
                      species='',
                      method='',
                      pdb='',
                      resolution='',
                      pref_chain='A',
                      state='Inactive',
                      deg_act='-',
                      family='-',
                      subtype='-',
                      function=''):
        cols = self.table.columns
        row = [uniprot, receptor_family, cl, species, method, pdb, resolution, pref_chain, state, deg_act, family,
               subtype, function]
        row_df = pd.DataFrame([row], columns=cols)
        return row_df

    # ======================================================================================================================

    def update_pdbs(self):
        updatepdbs(self.path_pdb)

    # ======================================================================================================================

def make_ref_file_csv(path=r'data/gpcrdb/residue_table.xlsx'):
    ref_table = pd.read_excel(path)
    csv_path = path.split('.')[0]+'.csv'
    ref_table.to_csv(csv_path, index=False)
    print("Wrote reference table data to {}!".format(csv_path))
    del ref_table
    
def make_ref_table(path=r'data/gpcrdb/residue_table.csv', to_csv=True):
    ref_table = pd.read_csv(path, index_col='GPCRdb(A)')
    indices = ref_table.index
    r_dict = {}
    for r, region in enumerate(regions):
        if region != 'H8':
            indices = list(ref_table.loc[region:regions[r+1]].index)[1:-1]
        else:
            indices = list(ref_table.loc[region:].index)[1:-1]
        r_dict.update({region: indices})
    ref_table = ref_table.drop(regions)
    columns = [x.replace('<', '').replace('i>', '').replace('/', '') 
               for x in list(ref_table.columns)]
    ref_table.columns = columns
    if to_csv:
        ref_table.to_csv(path.split('.')[0]+'.csv')
    return ref_table  
    
def get_ref_table(path=r'data/gpcrdb/residue_table.csv'):
    return pd.read_csv(path, index_col='GPCRdb(A)')

def query_ref_table_name(name, human=True):
    """
    This stuff is specific to the reference data table (human readable gene names > need to be parsed)
    """
    text = get_uniprot_gene_name_query(gene_name='+'.join(name), limit=8)
    split_list = []
    for t in text:
        text_ = (str(t.decode('utf-8')).split('\n'))
        split_list.append(text_)
    table = split_list[0]
    table_cols = table[0].split('\t')
    table_data = [x.split('\t') for x in table[1:-1]]
    table_df = pd.DataFrame(table_data, columns=table_cols)
    if human:
        table_df = table_df[table_df['Organism'].str.contains('Human')]
    seqs = split_list[1:]
    del text
    next_header = True
    seq_dict = {}
    info_dict = {}
    seq_iter = ''
    for s, seq in enumerate(seqs[0]):
        if len(seq) > 0:
            if ('>sp' in seq) & (s > 0):
                next_header = True
                seq_dict.update({uniprot_id: seq_iter})
                seq_iter = ''
            if next_header:
                header = seq.split('|')
                uniprot_id = header[1]
                info_dict.update({uniprot_id: header[1:]})
            if not next_header:
                seq_iter += seq
            next_header = False
    if human:
        valid_uids = list(table_df['Entry'])
        invalid_uids = [x for x in seq_dict.keys() if x not in valid_uids]
        for invalid_uid in invalid_uids:
            del seq_dict[invalid_uid]
            del info_dict[invalid_uid]
    return table_df, seq_dict, info_dict




def get_reference_annotation_dfs(ref_table):
    ref_dict = {}
    for r in range(len(ref_table.columns)):
        sample = ref_table.loc[:, columns[r]]
        sample = sample[sample!='-']
        grn = list(sample.to_dict().keys())
        res_resnr = list(sample.to_dict().values())
        res_type = [str(r)[0] for r in res_resnr]
        res_nr = [str(r)[1:] for r in res_resnr]
        assert len(res_type)==len(res_nr)==len(grn), print("number of entries for residue types, numbers and grns do not match!")
        df = pd.DataFrame([res_nr, res_type, grn]).T
        df.columns = ['res_nr', 'res_type', 'grn']
        ref_dict.update({columns[r]: df})
        del df
    return ref_dict

def save_grn_dicts(ref_dict, path='data/grn/ref/'):
    for ref in ref_dict.keys():
        df = ref_dict[ref]
        df.to_pickle(path+ref.replace(' ', '_')+'.pkl')
        
def get_grn_by_idx(ref_dict, idx=0):
    return ref_dict[list(ref_dict.keys())[idx]]

def get_grn_seq(ref_dict, idx=0):
    df = ref_dict[list(ref_dict.keys())[idx]]
    return df['res_type'].to_list()
    
def load_grn_dicts(path='data/grn/ref/'):
    """
    Loads all reference dictionaries containing generic residue numbers
    """
    files = get_filenames(path=path)
    grn_dict = {}
    for file in files:
        name = file[:-4]
        df = pd.read_pickle(path+file)
        grn_dict.update({name: df})
    return grn_dict