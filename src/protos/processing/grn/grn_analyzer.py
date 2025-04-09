import pandas as pd
from protos.processing.grn.grn_table_utils import *   # get_occs_of_cys, get_dists_from_cterminus, get_dists_from_h8
from protos.processing.aa_dicts import *


class GRNAnalyzer:
    def __init__(self,
                 grn_data: pd.DataFrame,
                 table_name = 'GRN default table',):
        self.table_name = table_name
        self.data = grn_data
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
        if on_reset_delete:
            self.features = pd.DataFrame(index=self.data.index.tolist())
            self.map = pd.DataFrame(columns=self.data.columns.tolist())
        self.data = self.load_grn_table()

    # =============================================================================================

    def get_entries(self, indices):
        pass

    def get_interval(self, interval):
        # check if interval is present in data
        return None

    def get_seq(self, idx):
        pass

    def get_seqs(self):
        pass

    def get_len(self, idx):
        pass

    def get_lens(self):
        return None

    def populate_map_features(self, mode='aminoacid', selection=['C']):
        if mode == 'aminoacid':
            print("applying amino acid map to GRN table, selected:", selection)
        if mode == 'characteristics':
            selection = AA_CHARACTERISTICS[selection]
            print("applying amino acid map to GRN table, selected:", selection)
        self.map = self.data.applymap(lambda x: 1 if x[0] in selection else 0)

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