from protos.io.fasta_utils import validate_fasta_format as validate_fasta_content, read_fasta, write_fasta
from protos.io.file_utils import get_filenames
from protos.loaders.uniprot_utils import *
from protos.io.paths.path_config import ProtosPaths, DataSource, join_path

import re
import os
from pathlib import Path
import pandas as pd
from tqdm import tqdm, trange
import itertools

COLS = ['uniprot', 'seq', 'gene', 'species', 'organism', 'info', 'dataset']


class UniprotDL:
    """
    This class downloads protein sequences from UniProt, processes and saves them to be used in alignments etc.
    
    The class follows the Protos path organization standards:
    - Sequences are stored in the standard fasta directory
    - Dataset files are stored in the metadata directory
    - The organization is based on file formats, not data sources
    """
    def __init__(self,
                 data_root=None,
                 dataset='gprot_2585',
                 reload=False,
                 limit=10,
                 # Legacy parameters for backward compatibility
                 path=None,
                 datasets='datasets/',
                 table='table/',
                 fastas='fastas/'):
        """
        Initialize the UniProt loader using standard path organization.
        
        Args:
            data_root: Root directory for data. If None, uses default from ProtosPaths.
            dataset: Name of the dataset file (without extension).
            reload: Whether to reload data even if it exists.
            limit: Maximum number of sequences to process.
            
            # Legacy parameters (deprecated, will be removed in future)
            path: Old-style base path (deprecated, use data_root instead)
            datasets: Old-style datasets subdirectory (deprecated)
            table: Old-style table subdirectory (deprecated)
            fastas: Old-style fastas subdirectory (deprecated)
        """
        # Handle legacy parameters for backward compatibility
        if path is not None:
            import warnings
            warnings.warn(
                "The 'path' parameter is deprecated and will be removed in a future version. "
                "Use 'data_root' instead.", 
                DeprecationWarning, 
                stacklevel=2
            )
            # If path is provided but not data_root, use path's parent as data_root
            if data_root is None:
                data_root = os.path.dirname(os.path.dirname(path))
        
        # Initialize ProtosPaths for standard path management
        self.paths = ProtosPaths(user_data_root=data_root, create_dirs=True)
        
        # Set dataset information
        self.dataset = dataset
        
        # Get standard directories from ProtosPaths
        self.sequence_dir = self.paths.get_processor_path('sequence', source=DataSource.USER)
        self.fasta_dir = self.paths.get_sequence_subdir_path('fasta_dir', source=DataSource.USER)
        self.metadata_dir = self.paths.get_sequence_subdir_path('metadata_dir', source=DataSource.USER)
        
        # Set specific paths for files
        self.dataset_file = join_path(self.metadata_dir, f"{dataset}.txt")
        self.gene_table_file = join_path(self.metadata_dir, f"gene_df_{dataset}.csv")
        
        # Create directories if they don't exist
        os.makedirs(self.metadata_dir, exist_ok=True)
        os.makedirs(self.fasta_dir, exist_ok=True)
        
        # Legacy path attributes for backward compatibility
        # These will be removed in a future version
        self._setup_legacy_paths(path, datasets, table, fastas)
        
        # Containers
        self.genes = []  # UniProt IDs from dataset
        self.data_list = []  # Intermediate storage
        self.gene_df = pd.DataFrame()  # Complete data
        
        # Settings
        self.reload = reload
        self.limit = limit
        
        # Load available dataset information
        self.available_datasets = self._get_available_datasets()
    
    def _setup_legacy_paths(self, path, datasets, table, fastas):
        """
        Set up legacy path attributes for backward compatibility.
        These attributes will be deprecated in a future version.
        
        Args:
            path: Base path for legacy structure
            datasets: Datasets subdirectory
            table: Table subdirectory
            fastas: Fastas subdirectory
        """
        if path is None:
            # Use data_root/uniprot as the legacy path
            path = join_path(self.paths.user_data_root, 'uniprot/')
        
        # Legacy path attributes
        self.path = path
        self.path_datasets = path + datasets
        self.path_dataset = self.path_datasets + self.dataset
        self.path_table = path + table
        self.path_fastas = path + fastas
        # Use metadata_dir/database for database files
        self.path_database = join_path(self.metadata_dir, 'database')
        
        # Create legacy directories for backward compatibility
        os.makedirs(self.path_datasets, exist_ok=True)
        os.makedirs(self.path_table, exist_ok=True)
        os.makedirs(self.path_fastas, exist_ok=True)
        os.makedirs(self.path_database, exist_ok=True)
        
        # Ensure backward compatibility by checking if dataset file exists
        if not os.path.exists(self.dataset_file) and os.path.exists(self.path_dataset + ".txt"):
            # Copy legacy dataset file to standard location
            import shutil
            shutil.copy2(self.path_dataset + ".txt", self.dataset_file)
    
    def _get_available_datasets(self):
        """
        Get list of available datasets in the metadata directory.
        
        Returns:
            List of dataset names (without extension)
        """
        if not os.path.exists(self.metadata_dir):
            return []
        
        # Look for .txt files in the metadata directory
        try:
            files = [f for f in os.listdir(self.metadata_dir) 
                    if f.endswith('.txt') and os.path.isfile(join_path(self.metadata_dir, f))]
            
            # Add legacy datasets for backward compatibility
            if os.path.exists(self.path_datasets):
                legacy_files = get_filenames(self.path_datasets)
                for f in legacy_files:
                    if f not in files:
                        files.append(f)
            
            # Remove extension to get dataset names
            return [os.path.splitext(f)[0] for f in files]
        except:
            # Fallback to legacy method
            return get_filenames(self.path_datasets)
    
    def load_dataset(self, ext='.txt'):
        """
        Load UniProt IDs from a dataset file in the standard metadata directory.
        
        Args:
            ext: File extension for the dataset file
            
        Returns:
            List of UniProt IDs
        """
        # Try to load from standard location first
        standard_file = self.dataset_file
        if os.path.exists(standard_file):
            with open(standard_file) as file:
                genes = [line.rstrip().split(' ') for line in file]
            self.genes = list(itertools.chain(*genes))
            return list(dict.fromkeys(self.genes))
        
        # Fallback to legacy location for backward compatibility
        legacy_file = self.path_dataset + ext
        if os.path.exists(legacy_file):
            with open(legacy_file) as file:
                genes = [line.rstrip().split(' ') for line in file]
            self.genes = list(itertools.chain(*genes))
            
            # Copy to standard location for future use
            os.makedirs(os.path.dirname(standard_file), exist_ok=True)
            import shutil
            shutil.copy2(legacy_file, standard_file)
            
            return list(dict.fromkeys(self.genes))
        
        raise FileNotFoundError(f"Dataset file not found at {standard_file} or {legacy_file}")
    
    def save_gene_table(self, ext='.csv', force_overwrite=True):
        """
        Save gene data table to the standard metadata directory.
        
        Args:
            ext: File extension for the table file
            force_overwrite: Whether to overwrite existing file
            
        Returns:
            Path to the saved file, or None if no data to save
        """
        if self.gene_df.empty:
            return None
        
        # Save to standard location
        os.makedirs(os.path.dirname(self.gene_table_file), exist_ok=True)
        self.gene_df.to_csv(self.gene_table_file, index=False)
        
        # Also save to legacy location for backward compatibility
        legacy_file = self.path_table + 'gene_df_' + self.dataset + ext
        os.makedirs(os.path.dirname(legacy_file), exist_ok=True)
        self.gene_df.to_csv(legacy_file, index=False)
        
        return self.gene_table_file
    
    def load_gene_table(self, dataset=None, ext='.csv', set_gene_df=True):
        """
        Load gene data table from the standard metadata directory.
        
        Args:
            dataset: Name of the dataset (default: use self.dataset)
            ext: File extension for the table file
            set_gene_df: Whether to set self.gene_df with loaded data
            
        Returns:
            Loaded DataFrame, or None if file not found
        """
        ds = dataset if dataset is not None else self.dataset
        
        # Try to load from standard location first
        standard_file = join_path(self.metadata_dir, f"gene_df_{ds}{ext}")
        if os.path.exists(standard_file):
            gene_df = pd.read_csv(standard_file)
            if set_gene_df:
                self.gene_df = gene_df
                return self.gene_df
            else:
                return gene_df
        
        # Fallback to legacy location for backward compatibility
        legacy_file = self.path_table + 'gene_df_' + ds + ext
        if os.path.exists(legacy_file):
            gene_df = pd.read_csv(legacy_file)
            
            # Copy to standard location for future use
            os.makedirs(os.path.dirname(standard_file), exist_ok=True)
            gene_df.to_csv(standard_file, index=False)
            
            if set_gene_df:
                self.gene_df = gene_df
                return self.gene_df
            else:
                return gene_df
        
        return None
    
    def _check_if_gene_df_exists(self, dataset, ext):
        """
        Check if gene DataFrame exists for the given dataset.
        
        Args:
            dataset: Name of the dataset
            ext: File extension for the table file
            
        Returns:
            True if file exists, False otherwise
        """
        ds = dataset if dataset is not None else self.dataset
        
        # Check standard location first
        standard_file = join_path(self.metadata_dir, f"gene_df_{ds}{ext}")
        if os.path.exists(standard_file):
            return True
        
        # Fallback to legacy location
        legacy_file = self.path_table + 'gene_df_' + ds + ext
        return os.path.exists(legacy_file)
    
    def save_uniprot_fasta(self, uniprot=None, mode='entry', ext='.fasta',
                          force_overwrite=True, validate=True, 
                          legacy_ext='.faa', dataset=None):
        """
        Save sequence data to FASTA files in the standard fasta directory.
        
        Args:
            uniprot: UniProt ID or list of IDs (default: use all in gene_df)
            mode: 'entry' for individual files, 'database' for combined file
            ext: File extension to use (should be .fasta for consistency)
            force_overwrite: Whether to overwrite existing files
            validate: Whether to validate FASTA format
            legacy_ext: Extension to use for legacy files (for backward compatibility)
            dataset: Name of the dataset (default: use self.dataset)
            
        Returns:
            List of saved file paths
        """
        ds = dataset if dataset is not None else self.dataset
        
        # Determine which IDs to save
        uids = list(self.gene_df['uniprot']) if uniprot is None else (
            uniprot if isinstance(uniprot, list) else [uniprot])
        
        # Verify IDs exist in dataframe
        if uids:
            missing_uids = [uid for uid in uids if uid not in list(self.gene_df['uniprot'])]
            if missing_uids:
                raise ValueError(f"UniProt IDs not found in data: {missing_uids}")
        
        saved_files = []
        
        if mode == 'entry':
            # Save individual FASTA files in the standard fasta directory
            for uid in uids:
                row = self.gene_df[self.gene_df['uniprot'] == uid]
                
                # Extract sequence and info
                seq = str(row['seq'].values[0])
                info = str(row['info'].values[0]) if 'info' in row.columns else f"UniProt:{uid}"
                
                # Create sequence dictionary
                sequence_dict = {uid: seq}
                
                # Determine file path in standard location
                standard_path = join_path(self.fasta_dir, f"{uid}{ext}")
                
                # Create directory if needed
                os.makedirs(os.path.dirname(standard_path), exist_ok=True)
                
                # Write using standard fasta_utils function
                write_fasta(sequence_dict, standard_path)
                saved_files.append(standard_path)
                
                # Also save to legacy location for backward compatibility
                legacy_path = join_path(self.path_fastas, f"{uid}{legacy_ext}")
                os.makedirs(os.path.dirname(legacy_path), exist_ok=True)
                
                # Use legacy format for backward compatibility
                fasta = self._parse_seq_to_fasta_format(info, seq)
                with open(legacy_path, 'w') as f:
                    for line in fasta:
                        f.write(line)
                        f.write('\n')
                
        else:  # database mode
            # Save a combined FASTA file in the metadata directory
            standard_path = join_path(self.metadata_dir, f"{ds}_db{ext}")
            
            # Create combined dictionary
            sequences = {}
            for uid in uids:
                row = self.gene_df[self.gene_df['uniprot'] == uid]
                seq = str(row['seq'].values[0])
                sequences[uid] = seq
            
            # Create directory if needed
            os.makedirs(os.path.dirname(standard_path), exist_ok=True)
            
            # Write using standard fasta_utils function
            write_fasta(sequences, standard_path)
            saved_files.append(standard_path)
            
            # Also save to legacy location for backward compatibility
            legacy_path = join_path(self.path_database, f"{ds}_db{legacy_ext}")
            os.makedirs(os.path.dirname(legacy_path), exist_ok=True)
            
            # Use legacy format for backward compatibility
            database = []
            for uid in uids:
                row = self.gene_df[self.gene_df['uniprot'] == uid]
                seq = str(row['seq'].values[0])
                info = str(row['info'].values[0]) if 'info' in row.columns else f"UniProt:{uid}"
                fasta = self._parse_seq_to_fasta_format(info, seq)
                database.extend(fasta)
            
            with open(legacy_path, 'w') as f:
                for line in database:
                    f.write(line)
                    f.write('\n')
        
        return saved_files
    
    def save_to_standard_location(self, uniprot_id, force_overwrite=True):
        """
        Save a UniProt sequence to the standard sequence location.
        
        Args:
            uniprot_id: UniProt ID to save
            force_overwrite: Whether to overwrite existing file
            
        Returns:
            Path to the saved file, or None if failed
        """
        # Verify the ID exists in the dataframe
        if uniprot_id not in list(self.gene_df['uniprot']):
            raise ValueError(f"UniProt ID not found: {uniprot_id}")
        
        # Save using the standard save_uniprot_fasta method
        saved_files = self.save_uniprot_fasta(uniprot=uniprot_id, mode='entry', force_overwrite=force_overwrite)
        
        # Return the path to the saved file
        return saved_files[0] if saved_files else None
    
    def _parse_seq_to_fasta_format(self, info, seq):
        """
        Create FASTA file format from sequence information.
        
        Args:
            info: Header information
            seq: Sequence string
            
        Returns:
            List of lines in FASTA format
        """
        seq_len = len(seq)
        n_lines = seq_len // 60
        lines = [seq[(i*60):min((i+1)*60, len(seq))] for i in range(n_lines + 1)]
        if len(lines[-1]) == 0:
            lines = lines[:-1]
        _info = '>'+info
        return [_info] + lines
    
    def download_gene_single_query(self, uniprot):
        """
        Downloads a list of properties from uniprot and safes them
        :param uniprot: Uniprot gene name
        :return:
        """
        data = get_uniprot(uniprot)
        
        # Extract sequence data from the DataFrame structure
        if data.empty:
            raise ValueError(f"No data returned for UniProt ID: {uniprot}")
            
        # Get sequence
        seq = data["Sequence"].iloc[0] if "Sequence" in data.columns else ""
        
        # Get gene name
        if "Entry Name" in data.columns:
            gene_name = data["Entry Name"].iloc[0]
        else:
            gene_name = f"{uniprot}_unknown"
            
        # Get organism information
        organism = data["Organism"].iloc[0] if "Organism" in data.columns else "Unknown"
        
        # Get protein name/description
        info = data["Protein names"].iloc[0] if "Protein names" in data.columns else f"Protein {uniprot}"
        
        # Parse species and gene from the gene name
        if "_" in gene_name:
            gene = gene_name.split('_')[0]
            species = gene_name.split('_')[1]
        else:
            gene = gene_name
            species = "unknown"
            
        return [uniprot, seq, gene_name, species, organism, info, self.dataset]

    def download_genes_single_query(self,
                                   batchsize=256,
                                   save=True,
                                   keep_in_mem=False,
                                   ext='.csv',
                                   force_overwrite=True):
        """
        Download and save multiple genes from UniProt.
        
        Args:
            batchsize: Number of genes to download in each batch
            save: Whether to save the gene table
            keep_in_mem: Whether to keep all data in memory
            ext: File extension for the gene table
            force_overwrite: Whether to overwrite existing data
            
        Returns:
            DataFrame with gene data
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
                    print(f"Could not download sequence of {gene} of dataset {self.dataset}.")
            gene_df_ = pd.DataFrame(self.data_list, columns=COLS)
            if (len(self.gene_df) == 0) or (keep_in_mem or save):
                self.gene_df = gene_df_
            else:
                self.gene_df.append(gene_df_, ignore_index=True)
            del gene_df_
            if save: # is always true otw there will be an error
                # saves or updates the table corresponding to the current dataset\!
                self.update_gene_df(ext='.csv', force_overwrite=True)
        
        return self.gene_df
    
    def update_gene_df(self, dataset=None, ext='.csv', force_overwrite=True):
        """
        Save gene DataFrame to both standard and legacy locations.
        
        Args:
            dataset: Name of the dataset (default: use self.dataset)
            ext: File extension for the table file
            force_overwrite: Whether to overwrite existing entries
            
        Returns:
            True if executed
        """
        ds = dataset if dataset is not None else self.dataset
        
        # Load existing gene_df from standard location if available
        standard_file = join_path(self.metadata_dir, f"gene_df_{ds}{ext}")
        if os.path.exists(standard_file):
            _gene_df = pd.read_csv(standard_file)
        else:
            # Try legacy location
            legacy_file = self.path_table + 'gene_df_' + ds + ext
            if os.path.exists(legacy_file):
                _gene_df = pd.read_csv(legacy_file)
            else:
                _gene_df = pd.DataFrame(columns=COLS)
        
        # Update the gene_df
        if force_overwrite:
            # Overwrite existing entries
            _gene_df = _gene_df[~_gene_df['uniprot'].isin(self.gene_df['uniprot'])]
        else:
            # Keep existing entries
            self.gene_df = self.gene_df[~self.gene_df['uniprot'].isin(_gene_df['uniprot'])]
        
        # Combine the datasets
        _gene_df = pd.concat([_gene_df, self.gene_df], ignore_index=True)
        
        # Save to standard location
        os.makedirs(os.path.dirname(standard_file), exist_ok=True)
        _gene_df.to_csv(standard_file, index=False)
        
        # Save to legacy location for backward compatibility
        legacy_file = self.path_table + 'gene_df_' + ds + ext
        os.makedirs(os.path.dirname(legacy_file), exist_ok=True)
        _gene_df.to_csv(legacy_file, index=False)
        
        return True

    def set_limit(self, limit):
        """
        Set a limit for the number of genes to process.
        
        Args:
            limit: Maximum number of genes to process
        """
        self.limit = limit

    def set_dataset(self, dataset='uniprot'):
        """
        Set the current dataset.
        
        Args:
            dataset: Name of the dataset
        """
        self.dataset = dataset
        self.dataset_file = join_path(self.metadata_dir, f"{dataset}.txt")
        self.gene_table_file = join_path(self.metadata_dir, f"gene_df_{dataset}.csv")
        self.path_dataset = self.path_datasets + self.dataset
    
    def remove_known_genes_from_genes(self, dataset=None, ext='.csv'):
        """
        Remove genes that are already in a dataset from the list of genes to download.
        
        Args:
            dataset: Name of the dataset (default: use self.dataset)
            ext: File extension for the table file
            
        Returns:
            True if genes were removed, False otherwise
        """
        ds = dataset if dataset is not None else self.dataset
        
        # Load gene table
        _gene_df = self.load_gene_df(dataset=ds, ext=ext, set_gene_df=False)
        
        if _gene_df is not None:
            # Get list of known genes
            known_genes = list(_gene_df['uniprot'])
            
            # Filter out known genes
            self.genes = [gene for gene in self.genes if gene not in known_genes]
            return True
        
        return False
    
    def get_gene(self, idx=0):
        """Get a gene from the gene DataFrame by index."""
        return self.gene_df.iloc[idx]

    def get_gene_seq(self, idx=0):
        """Get a gene sequence from the gene DataFrame by index."""
        return self.gene_df.iloc[idx]['seq']

    def __len__(self):
        """Return the number of genes in the gene DataFrame."""
        return len(self.gene_df)


# ======================================================================================================================

def load_grn_dicts(path='data/grn/ref/'):
    """
    Loads all reference dictionaries containing generic residue numbers. Note that the reference dictionary-generating
    functions are implemented under data_utils/gpcrdb/gpcrdb_loader.py \!
    """
    files = get_filenames(path=path)
    grn_dict = {}
    for file in files:
        name = file[:-4]
        df = pd.read_pickle(path + file)
        grn_dict.update({name: df})
    return grn_dict
