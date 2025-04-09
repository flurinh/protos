# should be able to use graphs
# should be able to define them in terms of grns
# should be able to ...

class LigandProcessor:
    def __init__(self, cif_processor):
        self.cif_processor = cif_processor

    def find_ligand_binder(self, ligand_name):
        has_ligand = []
        for pdb_id in self.cif_processor.pdb_ids:

            # Filter the data for HETATM records matching the specified ligand name
            ligand_data = self.cif_processor.data[(self.cif_processor.data['pdb_id'] == pdb_id) &
                                                  (self.cif_processor.data['group'] == 'HETATM') &
                                                  (self.cif_processor.data['res_name3l'] == ligand_name)]

            if not ligand_data.empty:
                has_ligand.append(pdb_id)
        return has_ligand

    def get_ligand(self, pdb_id, ligand_name):
        structure_data = self.cif_processor.get_structure_by_idx(pdb_id)
        ligand_data = structure_data[(structure_data['group'] == 'HETATM') &
                                     (structure_data['res_name3l'] == ligand_name)]
        return ligand_data

