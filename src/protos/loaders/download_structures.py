from Bio.PDB import PDBList
import os

def download_protein_structures(pdb_ids, target_folder='data/mmcif/'):
    pdbl = PDBList()
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    for pdb_id in pdb_ids:
        pdbl.retrieve_pdb_file(pdb_id, pdir=target_folder, file_format="mmCif")
