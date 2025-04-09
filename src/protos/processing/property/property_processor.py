import pandas as pd
import numpy as np
import os
import pickle
from sklearn.manifold import TSNE
import plotly.graph_objects as go
import plotly.colors as pc


class PropertyProcessor:
    def __init__(self, dataset=None, data_folder='data/properties/'):
        self.identity = pd.DataFrame()
        self.properties = pd.DataFrame()
        self.metadata = {}
        self.data_folder = data_folder
        self.dataset_name = dataset

        if dataset:
            self.load_dataset(dataset)
        else:
            self.available_identities = []
            self.available_properties = []

    def set_dataframe(self, dataframe_type, dataframe):
        """
        Sets the specified dataframe (identity or properties) and tracks its column datatypes.
        """
        self.check_protein_id_presence(dataframe)
        if dataframe_type == 'identity':
            self.identity = dataframe
            self.identity_dtypes = self._get_dtypes(dataframe)
            self.set_identities(self.identity.columns.tolist())
        elif dataframe_type == 'properties':
            self.properties = dataframe
            self.properties_dtypes = self._get_dtypes(dataframe)
            self.set_properties(self.properties.columns.tolist())
        else:
            raise ValueError("dataframe_type must be 'identity' or 'properties'.")

    def _get_dtypes(self, dataframe):
        """
        Returns a dictionary mapping column names to their datatypes in the given dataframe.
        """
        return dataframe.dtypes.apply(lambda x: x.name).to_dict()

    def _apply_filters(self, df, filters):
        """
        Apply filters to a dataframe.

        :param df: The DataFrame to filter.
        :param filters: Filtering conditions as a dictionary.
        :return: Filtered DataFrame.
        """
        for key, value in filters.items():
            if '__' in key:
                col, op = key.split('__', 1)
                if op == 'eq':
                    df = df[df[col] == value]
                elif op == 'gt':
                    df = df[df[col] > value]
                elif op == 'lt':
                    df = df[df[col] < value]
                elif op == 'contains':
                    df = df[df[col].astype(str).str.contains(value, na=False)]
                elif op == 'is_in':
                    df = df[df[col].isin(value)]
                elif op == 'not_na':
                    df = df[df[col].notna()]
            else:
                df = df[df[key] == value]
        return df

    def filter(self, df_type, filters, map_to_other=False, inplace=False):
        """
        General-purpose method to filter identity or property DataFrames and optionally map the filtering to the other DataFrame.

        :param df_type: 'identity' or 'properties' to indicate which DataFrame to filter.
        :param filters: Filtering conditions as a dictionary.
        :param map_to_other: If True, apply the filtering based on protein_id to the other DataFrame.
        :param inplace: If True, update the DataFrame stored in the class.
        :return: Filtered DataFrame.
        """
        df = self.identity if df_type == 'identity' else self.properties
        filtered_df = self._apply_filters(df, filters)

        if map_to_other:
            other_df = self.properties if df_type == 'identity' else self.identity
            mapped_ids = filtered_df['protein_id'].unique()
            mapped_df = other_df[other_df['protein_id'].isin(mapped_ids)]

            if inplace:
                self.set_dataframe(df_type, filtered_df)
                self.set_dataframe('identity' if df_type == 'properties' else 'properties', mapped_df)
            return filtered_df, mapped_df

        if inplace:
            self.set_dataframe(df_type, filtered_df)
        return filtered_df

    def add_new_column(self, new_data, data_type='properties', filter_missing=False):
        """
        Enhanced to handle datatypes correctly when adding new data and filling missing values.
        """
        if not isinstance(new_data, pd.DataFrame):
            new_data = pd.DataFrame(new_data)

        if 'protein_id' not in new_data.columns:
            raise ValueError("The new data must contain a 'protein_id' column.")

        target_df = self.identity if data_type == 'identity' else self.properties
        dtypes_dict = self.identity_dtypes if data_type == 'identity' else self.properties_dtypes

        # Ensure new data columns match target dataframe datatypes, particularly for new columns
        for col in new_data.columns:
            if col in dtypes_dict:
                new_data[col] = new_data[col].astype(dtypes_dict[col])
            else:
                # Track new column datatypes
                dtypes_dict[col] = new_data[col].dtype

        # Merge new data into the target DataFrame
        updated_df = pd.merge(target_df, new_data, on='protein_id', how='left')

        # Handle missing data according to column datatypes
        for col, dtype in dtypes_dict.items():
            if dtype == 'object':
                updated_df[col].fillna('', inplace=True)
            else:
                updated_df[col].fillna(np.nan, inplace=True)

        if filter_missing:
            updated_df.dropna(subset=['protein_id'], inplace=True)
            self.synchronize_dataframes(updated_df['protein_id'].unique())

        # Update the dataframe and its datatypes tracking
        if data_type == 'identity':
            self.identity = updated_df
            self.identity_dtypes = dtypes_dict
        else:
            self.properties = updated_df
            self.properties_dtypes = dtypes_dict

        self.set_properties(self.properties.columns.tolist())
        self.set_identities(self.identity.columns.tolist())

    def synchronize_dataframes(self, common_ids=None):
        """
        Synchronizes the identity and properties DataFrames based on common protein IDs.

        :param common_ids: Array of common protein IDs. If None, computes the common IDs.
        """
        if common_ids is None:
            identity_ids = set(self.identity['protein_id'])
            properties_ids = set(self.properties['protein_id'])
            common_ids = identity_ids.intersection(properties_ids)

        self.identity = self.identity[self.identity['protein_id'].isin(common_ids)]
        self.properties = self.properties[self.properties['protein_id'].isin(common_ids)]

    def check_protein_id_presence(self, dataframe):
        if 'protein_id' not in dataframe.columns:
            raise ValueError("The dataframe must contain a 'protein_id' column.")

    def check_protein_id_consistency(self):
        """
        Checks if the protein IDs in both identity and properties DataFrames are consistent, with no missing elements.

        :return: Boolean indicating consistency (True) or inconsistency (False).
        """
        self.check_protein_id_presence(self.identity)
        self.check_protein_id_presence(self.properties)
        identity_ids = set(self.identity['protein_id'])
        properties_ids = set(self.properties['protein_id'])
        return identity_ids == properties_ids

    def filter_by_identity(self, filters, map_to_properties=False, inplace=False):
        return self.filter('identity', filters, map_to_properties, inplace)

    def filter_by_property(self, filters, map_to_identity=False, inplace=False):
        return self.filter('properties', filters, map_to_identity, inplace)

    def get_protein_ids(self):
        return set(self.identity['protein_id'])

    def retain_columns(self, data_type, columns_to_keep):
        # Ensure 'protein_id' is always included in the columns to keep
        columns_to_keep = set(columns_to_keep + ['protein_id'])

        # Select the appropriate dataframe
        if data_type == 'identity':
            dataframe = self.identity
        elif data_type == 'properties':
            dataframe = self.properties
        else:
            raise ValueError("dataframe_type must be 'identity' or 'properties'.")

        # Filter the columns
        columns_to_remove = [col for col in dataframe.columns if col not in columns_to_keep]
        dataframe.drop(columns=columns_to_remove, inplace=True)

        # Update the class attribute with the modified dataframe
        if data_type == 'identity':
            self.identity = dataframe
            self.identity_dtypes = self._get_dtypes(dataframe)  # Update datatypes
        else:
            self.properties = dataframe
            self.properties_dtypes = self._get_dtypes(dataframe)  # Update datatypes

        self.set_properties(self.properties.columns.tolist())
        self.set_identities(self.identity.columns.tolist())

    def change_column_dtype(self, dataframe_type, column_name, new_dtype):
        """
        Attempts to change the datatype of a specified column in the given dataframe ('identity' or 'properties').

        :param dataframe_type: A string indicating the dataframe type ('identity' or 'properties').
        :param column_name: The name of the column whose datatype is to be changed.
        :param new_dtype: The target datatype to convert the column to.
        """
        if dataframe_type == 'identity':
            target_df = self.identity
        elif dataframe_type == 'properties':
            target_df = self.properties
        else:
            raise ValueError("dataframe_type must be 'identity' or 'properties'.")

        try:
            target_df[column_name] = target_df[column_name].astype(new_dtype)
            # Update datatypes tracking
            if dataframe_type == 'identity':
                self.identity_dtypes[column_name] = new_dtype
            else:
                self.properties_dtypes[column_name] = new_dtype
            print(f"Column '{column_name}' in '{dataframe_type}' successfully converted to {new_dtype}.")
        except Exception as e:
            print(f"Failed to convert column '{column_name}' in '{dataframe_type}' to {new_dtype}. Error: {e}")

    def set_properties(self, property_list):
        self.available_properties = property_list

    def set_identities(self, identity_list):
        self.available_identities = identity_list

    def set_dataset_name(self, dataset_name):
        self.dataset_name = dataset_name

    def add_rows(self, dataframe_type, new_data):
        if dataframe_type not in ['identity', 'properties']:
            raise ValueError("dataframe_type must be 'identity' or 'properties'.")

        target_df = self.identity if dataframe_type == 'identity' else self.properties
        dtypes_dict = self.identity_dtypes if dataframe_type == 'identity' else self.properties_dtypes

        # Check for 'protein_id' presence in new_data
        self.check_protein_id_presence(new_data)

        # Ensure new_data columns match the existing dataframe columns
        for col in new_data.columns:
            if col in dtypes_dict:
                new_data[col] = new_data[col].astype(dtypes_dict[col])
            else:
                raise ValueError(f"Column {col} not found in the existing dataframe columns.")

        # Append new data
        combined_df = pd.concat([target_df, new_data], ignore_index=True)

        # Fill missing data according to column data types
        for col, dtype in dtypes_dict.items():
            if dtype == 'object':
                combined_df[col].fillna('', inplace=True)
            elif 'float' in dtype or 'int' in dtype:
                combined_df[col].fillna(0, inplace=True)
            else:
                combined_df[col].fillna(combined_df[col].dtype.type(), inplace=True)

        # Update the dataframe and its datatypes tracking
        if dataframe_type == 'identity':
            self.identity = combined_df
            self.identity_dtypes = self._get_dtypes(combined_df)  # Update datatypes
        else:
            self.properties = combined_df
            self.properties_dtypes = self._get_dtypes(combined_df)  # Update datatypes

        # Optionally recompile metadata after adding new rows
        self.compile_metadata()

    def load_properties(self, identity_path, properties_path):
        identity = pd.read_csv(identity_path)
        self.check_protein_id_presence(identity)
        del identity

        properties = pd.read_csv(properties_path)
        self.check_protein_id_presence(properties)
        del properties

        self.set_properties(self.properties.columns.tolist())
        self.set_identities(self.identity.columns.tolist())

        self.identity_dtypes = self._get_dtypes(self.identity)
        self.properties_dtypes = self._get_dtypes(self.properties)

    def save_properties(self, identity_path, properties_path):
        self.identity.to_csv(identity_path, index=False)
        self.properties.to_csv(properties_path, index=False)

    def initialize_empty_dataset(self, dataset, folder='data/fasta/processed'):
        """
        Initializes the PropertyProcessor with data from a CSV file.

        :param csv_file: Path to the CSV file containing 'protein_id' and 'description' columns.
        """
        # Load the dataset
        csv_file = os.path.join(folder, f"{dataset}_ids.csv")
        identity_df = pd.read_csv(csv_file)

        # Check for required columns
        if 'protein_id' not in identity_df.columns or 'description' not in identity_df.columns:
            raise ValueError("CSV must contain 'protein_id' and 'description' columns.")

        # Set the identity DataFrame
        self.set_dataframe('identity', identity_df)

        # Create an empty properties DataFrame with only 'protein_id'
        self.properties = pd.DataFrame({'protein_id': identity_df['protein_id']})

        # Compile metadata for both DataFrames
        self.compile_metadata()

        # Optionally save the dataset
        self.save_dataset(dataset)

    def compile_metadata(self):
        """
        Compiles enhanced metadata for both identity and properties dataframes,
        including column names, data types, and various statistics based on data types.
        """

        def calculate_statistics(df):
            stats = {}
            for col, dtype in df.dtypes.items():
                col_stats = {'data_type': dtype.name}
                if dtype == 'float64':
                    col_stats['mean'] = df[col].mean()
                    col_stats['std'] = df[col].std()
                    col_stats['missing_values'] = df[col].isnull().sum()
                elif dtype == 'int64':
                    col_stats['value_counts'] = df[col].value_counts().to_dict()
                elif dtype == 'object':
                    # Assuming the column could contain string data or lists
                    if df[col].apply(lambda x: isinstance(x, list)).any():
                        col_stats['average_list_length'] = df[col].dropna().apply(len).mean()
                    else:
                        col_stats['unique_values'] = df[col].nunique()
                        col_stats['missing_values'] = df[col].isnull().sum()
                # Handling for additional specific data types could be added here
                stats[col] = col_stats
            return stats

        self.metadata['identity'] = calculate_statistics(self.identity)
        self.metadata['properties'] = calculate_statistics(self.properties)

    def save_dataset(self, dataset_name):
        """
        Saves the identity, properties, and metadata to CSV files (for DataFrames)
        and metadata to a pickle file.
        """

        self.set_dataset_name(dataset_name)
        if not os.path.exists(self.data_folder):
            os.makedirs(self.data_folder)

        identity_path = os.path.join(self.data_folder, f"{dataset_name}_identity.csv")
        properties_path = os.path.join(self.data_folder, f"{dataset_name}_properties.csv")
        metadata_path = os.path.join(self.data_folder, f"{dataset_name}_metadata.pkl")

        # Saving identity and properties DataFrames
        self.identity.to_csv(identity_path, index=False)
        self.properties.to_csv(properties_path, index=False)

        # Saving metadata using pickle
        with open(metadata_path, 'wb') as metafile:
            pickle.dump(self.metadata, metafile)

    def load_dataset(self, dataset_name):
        """
        Loads the identity from CSV file and optionally loads the properties if the file exists.
        If properties file does not exist, initializes the properties DataFrame with only the protein_id column from identity.
        """
        identity_path = os.path.join(self.data_folder, f"{dataset_name}_identity.csv")
        properties_path = os.path.join(self.data_folder, f"{dataset_name}_properties.csv")
        metadata_path = os.path.join(self.data_folder, f"{dataset_name}_metadata.pkl")

        # Loading identity DataFrame
        if os.path.exists(identity_path):
            self.identity = pd.read_csv(identity_path)
            self.check_protein_id_presence(self.identity)  # Ensure the 'protein_id' column exists
        else:
            raise FileNotFoundError(f"Identity file '{identity_path}' not found.")

        # Loading properties DataFrame, if it exists
        no_props = True
        if os.path.exists(properties_path):
            self.properties = pd.read_csv(properties_path)
            self.check_protein_id_presence(self.properties)  # Ensure the 'protein_id' column exists
        else:
            print(
                f"Properties file '{properties_path}' not found. Initializing properties with 'protein_id' from identity.")
            # Initialize properties DataFrame with just the 'protein_id' column from identity
            self.properties = pd.DataFrame(self.identity['protein_id'])

        # Loading metadata from pickle, if available
        if os.path.exists(metadata_path):
            with open(metadata_path, 'rb') as metafile:
                self.metadata = pickle.load(metafile)
        else:
            print("Metadata file not found. Metadata not loaded.")
        if no_props:
            self.compile_metadata()

        # Set the dataset name
        self.set_dataset_name(dataset_name)

        # Update available properties and identities
        self.set_properties(self.properties.columns.tolist())
        self.set_identities(self.identity.columns.tolist())

        # Update the dtypes tracking
        self.identity_dtypes = self._get_dtypes(self.identity)
        self.properties_dtypes = self._get_dtypes(self.properties)

    def available_datasets(self):
        """
        Lists the available datasets in the data folder.
        Only datasets consisting of both identity and properties are listed.
        """
        dataset_files = os.listdir(self.data_folder)
        datasets = []
        for file in dataset_files:
            if file.endswith('_identity.csv') and file.replace('_identity.csv', '_properties.csv') in dataset_files:
                datasets.append(file.replace('_identity.csv', ''))
        return datasets

    def get_protein_data(self, protein_id):
        """
        Retrieves the identity and property data for a given protein_id.

        :param protein_id: The unique identifier for the protein.
        :return: A tuple containing two dictionaries, one for identity and one for properties.
        """
        identity_data = self.identity[self.identity['protein_id'] == protein_id].to_dict(orient='records')
        property_data = self.properties[self.properties['protein_id'] == protein_id].to_dict(orient='records')

        # Return the first (and only) record from each list, or None if not found
        identity_dict = identity_data[0] if identity_data else None
        property_dict = property_data[0] if property_data else None

        return {'identity': identity_dict, 'properties': property_dict}

    def get_multiple_protein_data(self, protein_ids):
        """
        Retrieves the identity and property data for a list of protein_ids.

        :param protein_ids: A list of unique identifiers for the proteins.
        :return: A dictionary where keys are protein_ids and values are tuples containing identity and property dictionaries.
        """
        result = {}
        for protein_id in protein_ids:
            protein_data = self.get_protein_data(protein_id)
            result[protein_id] = protein_data

        return result

    @staticmethod
    def filter_help():
        """
        Provides detailed help and explanations on how to define filters for filtering identity and properties DataFrames.

        :return: A help string with instructions.
        """
        help_text = """
        Detailed Filter Definitions:
        - 'column__eq': Filters rows where the column's value equals the specified value. 
          Example: {'name__eq': 'ProteinA'} selects rows where the 'name' column is 'ProteinA'.

        - 'column__gt': Filters rows where the column's value is greater than the specified value.
          Example: {'length__gt': 30} selects rows where the 'length' column is greater than 30.

        - 'column__lt': Filters rows where the column's value is less than the specified value.
          Example: {'length__lt': 60} selects rows where the 'length' column is less than 60.

        - 'column__contains': Filters rows where the column's string value contains the specified substring.
          Example: {'function__contains': 'ase'} selects rows where the 'function' column contains 'ase'.

        - 'column__is_in': Filters rows where the column's value is in a provided list.
          Example: {'protein_id__is_in': [1, 2, 3]} selects rows where the 'protein_id' is either 1, 2, or 3.

        - 'column__not_na': Filters rows where the column's value is not NA (not missing).
          Example: {'location_not_na': True} selects rows where the 'location' column is not missing.

        Note:
        - The filter keys are case-sensitive and must exactly match the column names in the DataFrames.
        - Multiple filters can be combined, applying an 'AND' logic between them.
        - The '__contains' filter is only applicable to string-type columns.
        - The '__not_na' filter expects a boolean True as its value.

        Usage:
        After initializing the PropertyProcessor and loading or setting the identity and properties DataFrames,
        use the 'filter_by_identity' and 'filter_by_property' methods to apply these filters.
        """
        print(help_text)

    def __repr__(self):
        if self.metadata:
            formatted_metadata = "Dataset Metadata:\n\n"

            def format_dict(d, indent=0):
                result = ""
                for key, value in d.items():
                    if isinstance(value, dict):
                        result += "    " * indent + f"{key}:\n" + format_dict(value, indent + 1) + "\n"
                    else:
                        result += "    " * indent + f"{key}: {value}\n"
                return result

            formatted_metadata += format_dict(self.metadata)
        else:
            formatted_metadata = "No metadata available. Try loading a dataset or set dataframes and compile metadata."
        return formatted_metadata

    def load_and_merge_datasets(self, dataset_names):
        merged_identity = pd.DataFrame()
        merged_properties = pd.DataFrame()

        for dataset_name in dataset_names:
            identity_path = os.path.join(self.data_folder, f"{dataset_name}_identity.csv")
            properties_path = os.path.join(self.data_folder, f"{dataset_name}_properties.csv")

            if os.path.exists(identity_path) and os.path.exists(properties_path):
                identity_df = pd.read_csv(identity_path)
                properties_df = pd.read_csv(properties_path)

                # Add a column to track the original dataset
                identity_df['original_dataset'] = dataset_name
                properties_df['original_dataset'] = dataset_name

                merged_identity = pd.concat([merged_identity, identity_df], ignore_index=True)
                merged_properties = pd.concat([merged_properties, properties_df], ignore_index=True)
            else:
                print(f"Warning: Dataset '{dataset_name}' not found. Skipping.")

        if merged_identity.empty or merged_properties.empty:
            raise ValueError("No valid datasets were found to merge.")

        # Update the class attributes with the merged data
        self.identity = merged_identity
        self.properties = merged_properties

        # Update datatypes
        self.identity_dtypes = self._get_dtypes(self.identity)
        self.properties_dtypes = self._get_dtypes(self.properties)

        # Update available properties and identities
        self.set_properties(self.properties.columns.tolist())
        self.set_identities(self.identity.columns.tolist())

        # Compile metadata for the merged dataset
        self.compile_metadata()

        # Set a new dataset name for the merged dataset
        merged_dataset_name = "_".join(dataset_names)
        self.set_dataset_name(f"merged_{merged_dataset_name}")

        print(f"Merged {len(dataset_names)} datasets. New dataset name: {self.dataset_name}")


def tsne_visualization(embeddings, labels=None, emb_ids=None, shapes=None):
    # Convert embeddings to a list for TSNE processing
    emb_list = np.array([emb for emb in embeddings.values()])
    # Check if the embeddings are already flat
    if emb_list.ndim > 2:
        print("Warning: Embeddings have more than 2 dimensions and will be averaged for visualization.")
        emb_list = emb_list.mean(axis=1)  # Averaging if the embeddings are multi-dimensional

    # t-SNE Dimensionality Reduction
    tsne = TSNE(n_components=2, random_state=42)
    tsne_results = tsne.fit_transform(emb_list)

    # Prepare plot
    fig = go.Figure()

    # Define shape mappings
    shape_map = {0: 'circle', 1: 'star', 2: 'x', 3: 'triangle-up', 4: 'square'}

    if labels is not None:
        # Ensure labels list matches the number of embeddings
        if len(labels) != len(emb_list):
            raise ValueError("Length of labels does not match the number of embeddings.")

        # Convert labels to numeric values if they're not already
        numeric_labels = np.array(labels, dtype=float)

        # Create a continuous color scale
        colorscale = pc.sequential.Viridis

        # Normalize the labels to [0, 1] range for color mapping
        label_min, label_max = numeric_labels.min(), numeric_labels.max()
        normalized_labels = (numeric_labels - label_min) / (label_max - label_min)

        # Process shapes if provided
        if shapes is not None:
            if len(shapes) != len(emb_list):
                raise ValueError("Length of shapes does not match the number of embeddings.")
            unique_shapes = sorted(set(shapes))
            shape_indices = [unique_shapes.index(s) for s in shapes]

        # Create traces for each shape category
        for i, shape in enumerate(unique_shapes) if shapes is not None else [(-1, None)]:
            mask = np.array(shape_indices) == i if shapes is not None else np.ones(len(emb_list), dtype=bool)

            hovertexts = [f"ID: {emb_ids[j] if emb_ids else ''}, Value: {labels[j]}, Shape: {shape}"
                          for j in range(len(labels)) if mask[j]] if emb_ids else None

            fig.add_trace(go.Scatter(
                x=tsne_results[mask, 0],
                y=tsne_results[mask, 1],
                mode='markers',
                marker=dict(
                    color=normalized_labels[mask],
                    colorscale=colorscale,
                    colorbar=dict(title="Label Value"),
                    size=4,
                    symbol=shape_map.get(shape, 'circle') if shapes is not None else 'circle',
                    showscale=i == 0,  # Show colorbar only for the first trace
                ),
                text=hovertexts,
                hoverinfo='text' if emb_ids else 'none',
                name=f'Shape {shape}' if shapes is not None else 'All points'
            ))
    else:
        # If no labels are provided, plot all points in blue
        hovertexts = [emb_ids[i] if emb_ids else '' for i in range(len(emb_list))] if emb_ids else None
        fig.add_trace(go.Scatter(
            x=tsne_results[:, 0],
            y=tsne_results[:, 1],
            mode='markers',
            marker=dict(color='blue', size=8),
            text=hovertexts,
            hoverinfo='text' if emb_ids else 'none'
        ))

    # Update layout
    fig.update_layout(
        title='t-SNE visualization of embeddings',
        xaxis_title='t-SNE 1',
        yaxis_title='t-SNE 2',
        hovermode="closest"
    )
    fig.show()