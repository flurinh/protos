from protos.io.file_utils import *
from protos.loaders.uniprot_utils import *

import re
from pathlib import Path
import pandas as pd
from tqdm import tqdm, trange
import itertools

COLS = ['uniprot', 'seq', 'gene', 'species', 'organism', 'info', 'dataset']


class UniprotDL:
    """
    This class downloads gpcr & G protein sequences, processes and saves them to be used in alignments etc.
    ## https://github.com/protwis/gpcrdb_data/tree/master/protein_data
    """
    def __init__(self,
                 path='data/uniprot/',
                 datasets='datasets/',
                 dataset='gprot_2585',
                 table='table/',
                 fastas='fastas/',
                 reload=False,
                 limit=10):
        # Paths
        self.path = path
        self.path_datasets = self.path + datasets
        self.dataset = dataset
        self.path_dataset = self.path_datasets + self.dataset
        self.path_table = path + table
        self.path_fastas = path + fastas
        self.path_database = 'data/Blast/db/'

        # Containers
        self.genes = []  # contains the list of all genes in the datasets (defined by txt files with uniprot ids)
        self.data_list = []  # intermediate storage that contains information and gene sequences for batches
        self.gene_df = pd.DataFrame()  # complete

        # Information
        self.reload = reload
        self.limit = limit
        self.available_datasets = get_filenames(self.path_datasets)
        self.available_tables = get_filenames(self.path_table)

        # How do I store my data: 1) gene_list (dataset), 2) gene_df (full data), 3) individual fastas, 4) full fasta

        if not Path(self.path).exists():
            raise OSError(f"Directory {self.path} not found")

    #############
    #  SETTERS  #################################################################################################
    #############

    def set_limit(self, limit):
        """
        set a limit (number) of genes to process before interrupt
        """
        self.limit = limit

    def set_dataset(self, dataset='uniprot'):
        """
        Using this 
        """
        self.dataset = dataset
        self.path_dataset = self.path_datasets + self.dataset
        print("Using dataset(s): {} (path: {})".format(self.dataset, self.path_dataset))

    #############
    #  LOADERS  #################################################################################################
    #############

    def _check_if_gene_df_exists(self, dataset, ext):
        ds = dataset if dataset is not None else self.dataset
        print("dataset:", ds)
        if 'gene_df_' + ds + ext in self.available_tables:
            return True
        else:
            return False

    def load_gene_df(self, dataset=None, ext='.csv', set_gene_df=True):
        """
        This function loads existing gene dataframes for the current or specified dataset
        :param dataset: dataset to load the gene df for, if None use the default dataset (self.dataset)
        :param ext: what file extension is the table stored in?
        :param set_gene_df: whether the table should be safed as the working version or just return it
        :return: dataframe if it exists, None if it does not exist
        """
        ds = dataset if dataset is not None else self.dataset
        if self._check_if_gene_df_exists(ds, ext):
            # print("loading specified dataset [{}]".format(ds))
            gene_df = pd.read_csv(self.path_table + 'gene_df_' + ds + ext)
        else:
            # print("loading dataset [{}] failed (not available)".format(ds))
            return None
        if set_gene_df:
            self.gene_df = gene_df
            del gene_df
            return self.gene_df
        else:
            return gene_df

    def load_dataset(self, ext='.txt'):
        filename = self.path_dataset + ext
        with open(filename) as file:
            genes = [line.rstrip().split(' ') for line in file]
        self.genes = list(itertools.chain(*genes))
        return list(dict.fromkeys(self.genes))

    def remove_known_genes_from_genes(self, dataset=None, ext='.csv'):
        """
        This function removes uniprot ids from the list 'genes' that is used to download new uniprot genes.
        I.e. we exclude genes from our list that should not be downloaded - this is meant to work around 're-downloading'
        the same genes multiple times.
        :param dataset: optional name of a dataset to use as set of known genes
        :param ext: file extension
        :return:
        """
        ds = dataset if dataset is not None else self.dataset
        print("Removing genes known to the dataset [{}]!".format(ds))
        if self._check_if_gene_df_exists(dataset=ds, ext=ext):
            _gene_df = self.load_gene_df(dataset=None, ext='.csv', set_gene_df=False)
            known_genes = list(_gene_df['uniprot'])
            print("known (tobe deleted) genes: ", [y for y in self.genes if y in known_genes])
            self.genes = [x for x in self.genes if x not in known_genes]
            return True
        else:
            return False

    def _add_unknown_genes_to_genes(self, dataset=None, ext='.csv'):
        """
        This function adds uniprot ids form a gene_df that are NOT in the current working list of genes!
        Usecase: E.g. augmenting a dataset
        :param dataset: optional name of a dataset to use as set of known genes
        :param ext: file extension
        :return:
        """
        _gene_df = self.load_gene_df(dataset=None, ext='.csv', set_gene_df=False)
        known_genes = list(_gene_df['uniprot'])
        # add unkown genes
        self.genes = list(set(known_genes + self.genes))
        del _gene_df

    ################
    #  DOWNLADERS  ##############################################################################################
    ################

    def download_genes_single_query(self,
                                    batchsize=256,
                                    save=True,
                                    keep_in_mem=False,
                                    ext='.csv',
                                    force_overwrite=True):
        """
        :param batchsize: number of uniprot to be downloaded then saved before downloading continues - 
                          make sure breaking does not ruin the process
        :param save:      True if uniprot should be saved
        :param keep_in_mem:
                          True if we want to load the whole gene_df into memory
        """
        if self.limit is None:
            self.limit = len(self.genes)
        n_batches = self.limit // batchsize
        for b in trange(n_batches):
            for g, gene in enumerate(tqdm(self.genes[b * batchsize:(b + 1) * batchsize])):
                try:
                    print("getting new gene: ",gene)
                    self.data_list.append(self.download_gene_single_query(gene))
                except:
                    print("Could not download sequence of {} of dataset {}.".format(gene, self.dataset))
            gene_df_ = pd.DataFrame(self.data_list, columns=COLS)
            if (len(self.gene_df) == 0) or (keep_in_mem or save):
                self.gene_df = gene_df_
            else:
                self.gene_df.append(gene_df_, ignore_index=True)
            del gene_df_
            if save: # is always true otw there will be an error
                # saves or updates the table corresponding to the current dataset!
                self.update_gene_df(ext='.csv', force_overwrite=True)

    def download_gene_single_query(self, uniprot):
        """
        Downloads a list of properties from uniprot and safes them
        :param uniprot: Uniprot gene name
        :return:
        """
        data = get_uniprot(uniprot)
        seq = str(data.seq)
        gene_name = data.id[data.id.rfind('|') + 1:]
        info = data.description
        organism = info[info.index('OS') + 3:info.index('OX') - 1]
        species = gene_name.split('_')[1]
        gene = gene_name.split('_')[0]
        return [uniprot, seq, gene_name, species, organism, info, self.dataset]

    def download_genes_batch_query(self):
        """
        TODO: use a batch retrieval system: https://www.uniprot.org/help/api_queries // https://www.uniprot.org/help/api
        :return:
        """
        pass

    ############
    #  SAVERS  ####################################################################################
    ############

    def update_gene_df(self, dataset=None, ext='.csv', force_overwrite=True):
        """
        saves a gene dataframe (of a specific dataset)
        :param ext: defines file format (atm csv) -> should be really fast!
        :param force_overwrite: -> whether to update/not update existing uniprot entries
        :return: True if executed
        """
        ds = dataset if dataset is not None else self.dataset
        if ext == '.csv':
            if self._check_if_gene_df_exists(dataset=ds, ext=ext):
                _gene_df = self.load_gene_df(dataset=ds, ext=ext, set_gene_df=False)
            else:
                _gene_df = pd.DataFrame(columns=COLS)
            if force_overwrite:
                # Overwrite existing entries in stored gene_df with new data
                print("Gene df stored:", _gene_df)
                print("find copies:",list(_gene_df[_gene_df['uniprot'].isin(list(self.gene_df['uniprot']))].index))
                _gene_df.drop(index=list(_gene_df[_gene_df['uniprot'].isin(list(self.gene_df['uniprot']))].index),
                              inplace=True)
            else:
                # Do not overwrite existing entries!
                self.gene_df.drop(index=self.gene_df['uniprot'].isin(list(_gene_df['uniprot'])).index, inplace=True)
            print(len(_gene_df), len(self.gene_df))
            _gene_df = pd.concat([_gene_df, self.gene_df], ignore_index=True, axis=0)

            _gene_df.to_csv(self.path_table + 'gene_df_' + ds + ext, index=False)
            del _gene_df

    def save_uniprot_fasta(self, dataset=None, uniprot=None, mode='entry', ext='.faa',
                           force_overwrite=True, validate=True):
        """
        Save sequence data to fasta file
        :param dataset: which dataset
        :param uniprot:
        :param mode: 'entry' or 'database' -> either save all at once in a file or single fasta files
        :param ext:
        :return:
        """
        ds = dataset if dataset is not None else self.dataset
        if ds != self.dataset:
            # Load dataset to make fastas
            _gene_df = self.load_gene_df(dataset=ds, set_gene_df=False)
        else:
            # Use working memory dataset
            _gene_df = self.gene_df
        uids = list(_gene_df['uniprot']) if uniprot is None else uniprot if isinstance(uniprot, list) else [uniprot]
        assert len([uid for uid in uids if uid in list(_gene_df)]) == 0, \
            print("not all specified uniprot ids present in table!")
        if mode == 'entry':
            fasta_dict = {}
        else:
            database = []
        for uid in uids:
            row = _gene_df[_gene_df['uniprot']==uid]
            # This looks very ugly, somewhere I parsed a list to string and now am stuck with string literals
            info = str(row['info'].values).replace('\n', '')[2:-2]
            seq = str(row['seq'].values).replace(' ', '').replace('\n', '').replace('\'', '')\
                .replace('[','').replace(']','')
            fasta = self._parse_seq_to_fasta_format(info, seq)
            if validate:
                if not validate_fasta_format(fasta):
                    print("Found invalid Fasta format:", fasta)
                    print("\n\n\n")
            if mode == 'entry':
                fasta_dict.update({uid: fasta})
            else:
                database += fasta
        if mode == 'entry':
            for key in list(fasta_dict.keys()):
                filename = self.path_fastas + key + ext
                # TODO: check if file exists and comply with force_overwrite!
                with open(filename, 'w') as f:
                    for line in fasta_dict[key]:
                        f.write(line)
                        f.write('\n')
        else:
            filename = self.path_database + self.dataset + '_db' + ext
            with open(filename, 'w') as f:
                for line in database:
                    f.write(line)
                    f.write('\n')
        del _gene_df

    def _parse_seq_to_fasta_format(self, info, seq):
        """
        create fasta file format from
        :param info:
        :param seq:
        :return:
        """
        seq_len = len(seq)
        n_lines = seq_len // 60
        lines = [seq[(i*60):min((i+1)*60, len(seq))] for i in range(n_lines + 1)]
        if len(lines[-1]) == 0:
            lines = lines[:-1]
        _info = '>'+info
        return [_info] + lines







    #############
    #  GETTERS  #################################################################################################
    #############

    # Get data from Mem
    def get_gene(self, idx=0):
        return self.gene_df.iloc[idx]

    def get_gene_seq(self, idx=0):
        return self.gene_df.iloc[idx]['seq']

    def __len__(self):
        return len(self.gene_df)


# ======================================================================================================================

def load_grn_dicts(path='data/grn/ref/'):
    """
    Loads all reference dictionaries containing generic residue numbers. Note that the reference dictionary-generating
    functions are implemented under data_utils/gpcrdb/gpcrdb_loader.py !
    """
    files = get_filenames(path=path)
    grn_dict = {}
    for file in files:
        name = file[:-4]
        df = pd.read_pickle(path + file)
        grn_dict.update({name: df})
    return grn_dict
