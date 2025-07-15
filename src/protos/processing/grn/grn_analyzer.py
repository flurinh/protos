import pandas as pd
from typing import Optional, List, Dict, Union
from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.grn.grn_table_utils import *   # get_occs_of_cys, get_dists_from_cterminus, get_dists_from_h8
from protos.processing.grn.grn_utils import get_seq
from protos.processing.aa_dicts import *


class GRNAnalyzer:
    """Analyzer for GRN tables.
    
    This class provides analysis functionality for GRN tables loaded
    through a GRNProcessor instance.
    """
    
    def __init__(self,
                 grn_processor: Optional[GRNProcessor] = None,
                 grn_data: Optional[pd.DataFrame] = None,
                 table_name: str = 'GRN default table'):
        """Initialize GRN Analyzer.
        
        Args:
            grn_processor: GRNProcessor instance to use for data loading
            grn_data: Pre-loaded GRN data (if processor not provided)
            table_name: Name for the analysis
        """
        self.table_name = table_name
        self.processor = grn_processor
        
        # Use provided data or load from processor
        if grn_data is not None:
            self.data = grn_data
        elif self.processor is not None and self.processor.data is not None:
            self.data = self.processor.data
        else:
            raise ValueError("Either grn_data or grn_processor with loaded data must be provided")
            
        self.features = pd.DataFrame(index=self.data.index.tolist())
        self.map = pd.DataFrame(columns=self.data.columns.tolist())

    def __len__(self):
        return len(self.data)

    def __str__(self):
        return self.table_name + ' ' + str(len(self))

    # =============================================================================================

    def apply_interval(self, grn_interval):
        self.data = self.data[grn_interval]
        self.map = self.map[grn_interval]

    def apply_entry_selection(self, indices):
        self.data = None

    def reset_data(self, on_reset_delete=True):
        """Reset data from processor."""
        if self.processor is None:
            raise ValueError("Cannot reset data without a processor instance")
            
        if on_reset_delete:
            self.features = pd.DataFrame(index=self.data.index.tolist())
            self.map = pd.DataFrame(columns=self.data.columns.tolist())
            
        # Reload data from processor
        if self.processor.data is not None:
            self.data = self.processor.data
        else:
            raise ValueError("Processor has no loaded data")

    # =============================================================================================

    def get_entries(self, indices: List[str]) -> pd.DataFrame:
        """Get specific entries from the data."""
        return self.data.loc[indices]

    def get_interval(self, interval: List[str]) -> pd.DataFrame:
        """Get data for specific GRN interval.
        
        Args:
            interval: List of GRN positions to include
            
        Returns:
            DataFrame with only specified columns
        """
        # Filter to columns that exist in data
        existing_cols = [col for col in interval if col in self.data.columns]
        if not existing_cols:
            raise ValueError(f"None of the requested columns {interval} exist in data")
        return self.data[existing_cols]

    def get_seq(self, idx: str) -> str:
        """Get sequence for a specific index."""
        return get_seq(idx, self.data)

    def get_seqs(self) -> Dict[str, str]:
        """Get all sequences."""
        return {idx: get_seq(idx, self.data) for idx in self.data.index}

    def get_len(self, idx: str) -> int:
        """Get length of sequence for a specific index."""
        return len(get_seq(idx, self.data))

    def get_lens(self) -> pd.Series:
        """Get lengths of all sequences."""
        return pd.Series({idx: self.get_len(idx) for idx in self.data.index})

    def populate_map_features(self, mode='aminoacid', selection=['C']):
        if mode == 'aminoacid':
            print("applying amino acid map to GRN table, selected:", selection)
        if mode == 'characteristics':
            selection = AA_CHARACTERISTICS[selection]
            print("applying amino acid map to GRN table, selected:", selection)
        self.map = self.data.map(lambda x: 1 if x[0] in selection else 0)

    def populate_length_features(self, grn_interval: list = []):
        lens = []
        for uid in self.data.index.tolist():
            if len(grn_interval) > 0:
                # this func is terribly slow... we should use apply/map to figure this out
                self.apply_interval(grn_interval)
                lens.append(len(get_seq(uid, self.data[grn_interval])))
            else:
                lens.append(len(get_seq(uid, self.data)))
        if len(grn_interval) > 0:
            self.features['length_' + grn_interval[0] + '_' + grn_interval[-1]] = lens
        else:
            self.features['length'] = lens
        return self.features

    def populate_dist_features(self, origin_grn: str = '8x52'):
        # this is applied to a binary map
        # we need mul
        # the result is a list for each index containing the distance at which the map is 1
        pass

    def find_charged_patches(self, bin_size, n_min, grn_interval: list = []):
        if len(grn_interval) > 0:
            self.apply_interval(grn_interval)
        pass

    def vectorize_seq(self, model, idx):
        seq = self.get_seq(idx)
        vector = model(seq)
        return vector

    def vectorize_seqs(self, model, feature_name = 'embedding'):
        while(feature_name in self.features.columns.tolist()):
            if '_' in feature_name:
                name, colidx = feature_name.split('_')
                feature_name = name + '_' + str(int(colidx) + 1)
            else:
                feature_name = feature_name + '_0'
        self.features[feature_name] = None  # apply the model to the sequence and vectorize it

    def aggr_map_features(self, method='sum'):
        if method == 'sum':
            return self.map.sum()
        if method == 'mean':
            return self.map.mean()
        if method == 'norm':
            return self.map.mean() / self.map.std()