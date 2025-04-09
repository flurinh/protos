import pandas as pd
from typing import Optional

from protos.processing.property.aff_processor import AffinityProcessor
from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.property.gene_processor import GeneProcessor
from protos.processing.structure.struct_base_processor import CifBaseProcessor


class MappingProcessor:
    def __init__(self,
                 grnp: Optional[GRNProcessor] = None,
                 ap: Optional[AffinityProcessor] = None,
                 gp: Optional[GeneProcessor] = None,
                 sp: Optional[CifBaseProcessor] = None,
                 populate_maps=True
                 ):
        self.data = {
            'grnp': grnp.data if grnp else pd.DataFrame(),
            'ap': ap.data if ap else pd.DataFrame(),
            'gp': gp.data if gp else pd.DataFrame(),
            'sp': sp.data if sp else pd.DataFrame(),
        }
        self.populate_maps = populate_maps
        self.maps = {}
        mapping_cols = ['query_name', 'grn_entry', 'uniprot', 'uni_entry_name', 'genes', 'organism', 'len']
        self.main_map = pd.DataFrame(columns=mapping_cols)

    def _merge_dataframes(self, df1, df2):
        # assuming 'id' is the common column
        merged_df = pd.merge(df1, df2, on='id')
        return merged_df

    def _map_grn_to_gp(self):
        return self._merge_dataframes(self.data['grnp'], self.data['gp'])

    def _map_ap_to_gp(self):
        return self._merge_dataframes(self.data['ap'], self.data['gp'])

    def _map_sp_to_gp(self):
        return self._merge_dataframes(self.data['sp'], self.data['gp'])

    def _map_grn_to_sp(self):
        return self._merge_dataframes(self.data['grnp'], self.data['sp'])

    def _map_grn_to_ap(self):
        return self._merge_dataframes(self.data['grnp'], self.data['ap'])

    def _map_ap_to_sp(self):
        return self._merge_dataframes(self.data['ap'], self.data['sp'])

    def create_mapping(self, a='grnp', b='ap'):
        if ((a == 'gp') and (b == 'grnp')) or ((a == 'grnp') and (b == 'gp')):
            map_instance = self._map_grn_to_gp()
        elif ((a == 'gp') and (b == 'sp')) or ((a == 'sp') and (b == 'gp')):
            map_instance = self._map_sp_to_gp()
        elif ((a == 'gp') and (b == 'ap')) or ((a == 'ap') and (b == 'gp')):
            map_instance = self._map_ap_to_gp()
        elif ((a == 'grnp') and (b == 'sp')) or ((a == 'sp') and (b == 'grnp')):
            map_instance = self._map_grn_to_sp()
        elif ((a == 'grnp') and (b == 'ap')) or ((a == 'ap') and (b == 'grnp')):
            map_instance = self._map_grn_to_ap()
        elif ((a == 'ap') and (b == 'sp')) or ((a == 'sp') and (b == 'ap')):
            map_instance = self._map_ap_to_sp()
        else:
            raise ValueError("Invalid arguments")
        return map_instance

