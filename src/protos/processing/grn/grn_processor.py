from protos.processing.grn.grn_utils import *
import os
import pandas as pd
import plotly.graph_objs as go


class GRNProcessor:
    def __init__(self,
                 dataset=None,
                 path='data/grn/ref/',
                 preload=True):
        if dataset is None:
            print("dataset cannot be None. Using ref")
            dataset = 'ref'
            # Keep the user-provided path instead of overwriting it

        self.path = path
        self.dataset = dataset if isinstance(dataset, str) else 'merged_'+'_'.join(dataset)
        self.filename = self.dataset + '.csv'
        self.ids = []
        self.grns = []
        self.features = pd.DataFrame()
        self.maps = {}

        if preload:
            if isinstance(dataset, list):
                self.load_and_merge_grn_tables(datasets=dataset, path=path)
            else:
                self.load_grn_table(dataset=self.dataset)
            self.map = pd.DataFrame(columns=self.grns)
        else:
            self.data = pd.DataFrame()
            self.map = pd.DataFrame()

    def list_available_datasets(self):
        files = os.listdir(self.path)  # List all files in the specified path
        datasets = [f.split('.')[0] for f in files if f.endswith('.csv')]  # Filter for .csv files and remove extension
        return datasets

    def save_grn_table(self, dataset=None):
        if dataset is not None:
            self.dataset = dataset
            self.filename = dataset + '.csv'
        elif self.dataset is None:
            raise ValueError("No dataset specified.")
        file_path = os.path.join(self.path, self.dataset + '.csv')
        self.data.to_csv(file_path, index=True)

    def load_grn_table(self, dataset=None, low_memory=False, remove_duplicates=True, path=None):
        if path is not None:
            self.path = path
        if dataset is not None:
            self.dataset = dataset
            self.filename = dataset + '.csv'
        file_path = os.path.join(self.path, self.filename)
        print('Loading data from', file_path)
        df = pd.read_csv(file_path, low_memory=low_memory)
        df = df.rename(columns={'Unnamed: 0': 'uen'})
        df = df.sort_values(by='uen', ascending=True)
        df = df.set_index('uen').fillna('-')
        self.data = df
        
        if remove_duplicates:
            self.remove_duplicate_ids()
        self.ids = self.data.index.tolist()
        self.grns = self.data.columns.tolist()
        self.grns = sort_grns_str(self.grns)
        self.data = self.data[self.grns]
        self.data.fillna('-', inplace=True)
        return self.data

    def get_available_grn_tables(self):
        files = os.listdir(self.path)
        return [f.split('.')[0] for f in files if f.endswith('.csv')]

    def get_seq_dict(self):
        grns = list(self.data.columns)
        sorted_grns = sort_grns_str(grns)
        grn_table = self.data.loc[:, sorted_grns]
        seqs = [get_seq(idx, grn_table) for idx in self.ids]
        return dict(zip(self.ids, seqs))

    def remove_duplicate_ids(self):
        print("Removing duplicate IDs...")
        seen = set()
        duplics = [x for x in self.ids if x in seen or seen.add(x)]
        duplics = list(set(duplics))
        if len(duplics) > 0:
            print("Found {} duplicate IDs, removing. {}.".format(len(duplics), duplics))
            singles = [x for x in self.ids if x not in duplics]
            df1 = self.data.loc[singles, :]
            
            # Fix for the TypeError: Create a list of Series for concatenation
            df2_rows = []
            for dupli in duplics:
                try:
                    # Get the first occurrence of each duplicate
                    row = self.data.loc[dupli, :].iloc[0]
                    df2_rows.append(row)
                except (IndexError, AttributeError):
                    # Handle case where the duplicate ID doesn't exist in the data
                    continue
            
            # Only concatenate if there are rows to add
            if df2_rows:
                df2 = pd.DataFrame(df2_rows, index=[dupli for dupli in duplics if dupli in self.data.index])
                self.data = pd.concat([df1, df2])
            else:
                self.data = df1
                
            self.ids = self.data.index.tolist()

    def load_and_merge_grn_tables(self, datasets, path='data/grn/ref/', remove_duplicates=True):
        self.data = pd.DataFrame()
        tables = []
        for dataset in datasets:
            grn_table = self.load_grn_table(dataset=dataset, path=path, remove_duplicates=remove_duplicates)
            tables.append(grn_table)
        self.data = pd.concat(tables, axis=0)
        self.grns = sort_grns_str(self.data.columns.tolist())
        self.data = self.data[self.grns]
        self.data.fillna('-', inplace=True)
        self.ids = self.data.index.tolist()
        self.dataset = 'merged_'+'_'.join(datasets)
        return self.data

    def reset_data(self, reset_maps=False, reset_features=False):
        self.load_grn_table()
        if reset_features:
            self.features = pd.DataFrame(index=self.data.index.tolist())
        if reset_maps:
            self.maps = {}
        self.ids = self.data.index.tolist()
        self.grns = self.data.columns.tolist()

    def apply_interval(self, grn_interval, apply_to_maps=True):
        grn_interval = [col for col in grn_interval if col in self.grns]
        self.data = self.data[grn_interval]
        if apply_to_maps:
            self.apply_interval_to_map(grn_interval)
        self.grns = self.data.columns.tolist()

    def get_grn_dict(self, reset_data=False):
        if reset_data:
            self.reset_data()
        grns = list(self.data.columns)
        sorted_grns = sort_grns_str(grns)
        grn_table = self.data.loc[:, sorted_grns]
        residue_mask = grn_table.replace('-', pd.NA).notna()
        grn_dict = {uen: residue_mask.columns[residue_mask.loc[uen]].tolist() for uen in residue_mask.index}
        return grn_dict

    def filter_data_by_occurances(self, threshold):
        # Calculate the number of '-' entries in each column of the 'data' DataFrame
        non_existent_counts = (self.data == '-').sum()
        # Calculate the number of genes in each column
        gene_counts = len(self.data) - non_existent_counts
        # Get the column names where the gene count is greater than or equal to the threshold
        filtered_columns = gene_counts[gene_counts >= threshold].index
        # Update the 'data' DataFrame, the maps, and the 'grns' list using the 'apply_interval' function
        self.apply_interval(filtered_columns)

    def filter_by_ids(self, ids):
        ids = [idx for idx in ids if idx in self.ids]
        self.data = self.data.loc[ids, :]
        self.ids = ids

    def sort_columns(self):
        cols_unsorted = [parse_grn_str2float(x) for x in self.data.columns.tolist()]
        cols_sorted = [parse_grn_float2str(x) for x in sort_grns(cols_unsorted)]
        self.data = self.data.loc[:, cols_sorted]
        self.grns = self.data.columns.tolist()

    # Maps
    def get_maps(self):
        return list(self.maps.keys())

    def apply_interval_to_map(self, grn_interval):
        maps = self.get_maps()
        for map_key in maps:
            self.maps[map_key] = self.maps[map_key][grn_interval]

    def populate_map_features(self, grn_interval: list = [], map_name=None, mode='aminoacid', amino_acids=['C'],
                              window_size=5, threshold=2):
        if map_name is None:
            map_name = mode

        if map_name in self.maps:
            raise ValueError(f"A map with the name '{map_name}' already exists.")

        if len(grn_interval) > 0:
            self.apply_interval(grn_interval)
        if len(grn_interval) > 0:
            self.apply_interval(grn_interval)
        if mode == 'aminoacid':
            print("Applying amino acid map to GRN table, selected:", amino_acids)
            map = self.data.applymap(lambda x: 1 if x[0] in amino_acids else 0)
        elif mode == 'charged_patch':
            amino_acids = ['D', 'E', 'K', 'R', 'H']
            print("Finding charged patches in GRN table, selected:", amino_acids)
            print("Note that we use the full available sequence to find charged patches, incomplete sequences cannot"
                  "be noticed.")
            map_cols = list(self.data.columns)
            map = pd.DataFrame(columns=map_cols)
            seq_dict = self.get_seq_dict()
            for idx in self.ids:
                seq = seq_dict[idx]
                seq_map = [0] * len(seq)

                for i in range(len(seq) - window_size + 1):
                    window = seq[i:i + window_size]
                    charged_count = sum(aa in window for aa in amino_acids)
                    if charged_count >= threshold:
                        seq_map[i:i + window_size] = [1] * window_size

                # Update the map with the new seq_map
                seq_id = 0
                for _, grn in enumerate(map.columns):
                    if self.data.at[idx, grn] != '-':
                        map.at[idx, grn] = seq_map[seq_id]
                        seq_id += 1
        else:
            print('Mode {} not specified\! Please select from:\n{}'.format(
                mode, ['aminoacid', 'charged_patch']
            ))

        self.maps[map_name] = map

    def aggr_map_features(self, map_name, method='sum'):
        if map_name not in self.maps:
            raise ValueError(f"No map found with the name '{map_name}'.")

        map_data = self.maps[map_name]

        if method == 'sum':
            return map_data.sum(skipna=True)
        if method == 'mean':
            return map_data.mean(skipna=True)

    def populate_uen_features(self, reset_data=True):
        if reset_data:
            self.reset_data(reset_features=False, reset_maps=False)
        self.features['uen'] = self.features.index.tolist()

    def populate_length_features(self, grn_interval: list = [], name='', reset_data=True):
        if reset_data:
            self.reset_data(reset_features=False, reset_maps=False)
        lens = []
        if len(grn_interval) > 0:
            self.apply_interval(grn_interval, apply_to_maps=False)
        for idx in self.ids:
            lens.append(len(get_seq(idx, self.data)))
        col_name = 'length'
        if len(name) > 0:
            col_name = col_name + '_' + name
        self.features[col_name] = lens
        self.features.index = self.ids

    def get_occs_of_aa(self, grn_start='8x52', aa='C'):
        columns = self.data.columns.tolist()
        idx = columns.index(grn_start)
        cols = columns[idx:]
        results = {}
        indices = self.data.index.tolist()
        for i, idx in enumerate(indices):
            occs = 0
            for g, grn in enumerate(cols):
                x = self.data.loc[idx, grn]
                if x != '-':
                    if aa in x:
                        occs += 1
            results.update({idx: occs})
        return results

    def get_dists_to_aa(self, grn_start='8x52', aa='C'):
        columns = self.data.columns.tolist()
        idx = columns.index(grn_start)
        cols = columns[idx:]
        results = {}
        indices = self.data.index.tolist()
        for i, idx in enumerate(indices):
            res = []
            start = 0
            for g, grn in enumerate(cols):
                x = self.data.loc[idx, grn]
                if x != '-':
                    if g == 0:
                        start = int(x[1:])
                    if aa in x:
                        res.append(int(x[1:]) - start + 1)
            if len(res) >= 1:
                results.update({idx: res})
            else:
                results.update({idx: [0]})
        return results

    def apply_custom_function(self, custom_function, **kwargs):
        results = custom_function(self.data, **kwargs)
        self.update_features(results)

    def update_features(self, results):
        for col_name, values in results.items():
            self.features[col_name] = self.features.index.map(values)

    # Visualization
    def plot_map_features(self, map_name, method='sum', title='Aggregated Map Features', x_label='GRN',
                          show_gene_counts=False, color_regions=False):
        if map_name not in self.maps:
            raise ValueError(f"No map found with the name '{map_name}'.")

        if method == 'sum':
            aggr_map_data = self.aggr_map_features(map_name, 'sum')
            y_label = 'Sum'
            y_range = [0, len(aggr_map_data)]
        elif method == 'mean':
            aggr_map_data = self.aggr_map_features(map_name, 'mean')
            y_label = 'Mean'
            y_range = [0, 1]
        else:
            raise ValueError("Invalid aggregation method. Choose from 'sum' or 'mean'")

        colors = []
        if color_regions:
            for col in aggr_map_data.index:
                colors.append(map_grn_to_color(col))

        data = []
        for i in range(len(aggr_map_data) - 1):
            trace = go.Scatter(x=aggr_map_data.index[i:i + 2], y=aggr_map_data[i:i + 2], mode='lines',
                               line=dict(color=colors[i] if color_regions else 'black'),
                               showlegend=False)
            data.append(trace)

        layout = go.Layout(title=title,
                           xaxis=dict(title=x_label, showgrid=False, zeroline=False),
                           yaxis=dict(title=y_label, showgrid=False, zeroline=False, range=y_range),
                           yaxis2=dict(title='Gene Counts', overlaying='y', side='right', showgrid=False,
                                       zeroline=False),
                           legend=dict(orientation='h', x=0.5, xanchor='center', y=-0.3, yanchor='bottom'),
                           plot_bgcolor='white')

        if show_gene_counts:
            non_existent_counts = (self.data == '-').sum()
            gene_counts = len(self.data) - non_existent_counts
            trace2 = go.Scatter(x=gene_counts.index, y=gene_counts, mode='lines', name='Gene Counts', yaxis='y2',
                                line=dict(color='gray', dash='dot'))
            data.append(trace2)

        fig = go.Figure(data=data, layout=layout)
        fig.show()

    def generate_grn_mutant(self, idx, mutations):
        new_row = self.data.loc[idx, :].copy()  # Make sure to create a copy of the WT row

        # Process each mutation in the list
        for mutation in mutations:
            wt_aa = mutation[0]  # Extract wild-type amino acid
            mut_aa = mutation[-1]  # Extract mutant amino acid
            residue_num = mutation[1:-1]  # Extract residue number

            # Identify and update the column corresponding to the mutation
            for col, value in new_row.iteritems():
                if str(value).startswith(wt_aa + residue_num):  # Check if cell starts with <wt_aa><residue_num>
                    new_row[col] = mut_aa + residue_num  # Update to <mut_aa><residue_num>
                    break  # Assuming only one column per residue number, we can break after finding the match

        # Create a DataFrame from the new row for compatibility with appending to grnp.data
        row_df = pd.DataFrame([new_row])

        return row_df
