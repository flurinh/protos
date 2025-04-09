"""
def _map_grns_from_chain_to_structure(self, grn_dict):
    print("Starting GRN mapping...", list(grn_dict.keys()))  # Debugging statement

    # Iterate over each chain and its GRN mappings
    for pdb_chain, grn_mappings in grn_dict.items():
        pdb_id, chain_id = pdb_chain.split('_')
        print(f"Processing {pdb_chain}...")  # Debugging statement

        if pdb_chain in self.chain_dict_ca_atom_ids:
            ca_atom_ids = self.chain_dict_ca_atom_ids[pdb_chain]
            if pdb_chain in self.chain_dict:
                sequence = self.chain_dict[pdb_chain]

                for res_index, ca_atom_id in enumerate(ca_atom_ids):
                    if res_index < len(sequence):
                        gen_seq_ids = self.data.loc[(self.data['pdb_id'] == pdb_id) &
                                                    (self.data['auth_chain_id'] == chain_id) &
                                                    (self.data['atom_id'] == ca_atom_id), 'gen_seq_id']

                        if len(gen_seq_ids) == 0:
                            print(f"No gen_seq_id found for {pdb_id}, {chain_id}, {ca_atom_id}")  # Debug statement
                            continue

                        for gen_seq_id in gen_seq_ids:
                            residue = sequence[res_index]
                            residue_number = residue + str(int(res_index + 1))
                            if residue_number in grn_mappings:
                                grn = grn_mappings[residue_number]
                                self.data.loc[(self.data['pdb_id'] == pdb_id) &
                                              (self.data['auth_chain_id'] == chain_id) &
                                              (self.data['gen_seq_id'] == gen_seq_id), 'grn'] = grn
                            else:
                                print(f"Residue {residue_number} not in GRN mappings for {pdb_id}, "
                                      f"{chain_id}")

    print("GRN mapping completed.")

def assign_grns(self, grnp: GRNProcessor, aligner, min_score=0.4, load_alignment=True, verbose=1):
    self.data['seq_nr'] = ''
    self.data['grn'] = ''  # Initialize GRN column

    query_dict = self.chain_dict  # Prepare query sequences for alignment

    # Perform global alignment using mmseqs2
    ref_dict = grnp.get_seq_dict()
    alignment_df = mmseqs2_align2(query_dict, ref_dict)

    grn_dict = {}
    print(f"Using preprocessed alignments: {load_alignment}.")
    for pdb_chain_id, seq in self.chain_dict.items():
        print(f"Assigning GRNs to {pdb_chain_id}...")
        alignment_filename = f'{self.path_alignment_dir}/{pdb_chain_id}.txt'
        if load_alignment & os.path.isfile(alignment_filename):
            with open(alignment_filename, 'r') as f:
                grn_dict[pdb_chain_id] = json.load(f)
        else:
            alignment_df_result = alignment_df[(alignment_df['query_id'] == pdb_chain_id) &
                                               (alignment_df['sequence_identity'] >= min_score)]
            # Find the best alignment using the alm_df and filter for cutoff
            if len(alignment_df_result) == 0:
                print(f'Alignment insufficient {pdb_chain_id}!')
                continue

            # Get the best sequence alignment using blosum62
            try:
                best_match = alignment_df_result.iloc[0]['target_id']
                if verbose >= 1:
                    print("best_match", best_match)
                query_chain_dict = {pdb_chain_id: seq}
                ref_dict = {best_match: get_seq(best_match, grnp.data)}
                msa = msa_blosum62(query_chain_dict, ref_dict, aligner)
                alignment = msa[pdb_chain_id][2]
                row = grnp.data.loc[best_match]
                grn_list, rn_list, missing = expand_annotation(row, seq, alignment)
                grn_dict[pdb_chain_id] = dict(zip(rn_list, grn_list))

                # Save the grn_dict to a file
                with open(alignment_filename, 'w') as f:
                    json.dump(grn_dict[pdb_chain_id], f)
            except:
                print("Error producing the grn annotation of", pdb_chain_id)
    print(f"Mapping GRNs to the structure {pdb_chain_id}.")
    if len(grn_dict) > 0:
        self._map_grns_from_chain_to_structure(grn_dict)"""






"""def calculate_missing_ctail_grns(aligned_grns, missing_gene_numbers, query_gene_len, grns_float):
    ending_h8_float = grns_float[-1]
    print("LAST 10 std grns", grns_float[-10:])
    print("missing gene numbers", missing_gene_numbers)
    last_in_sequence = missing_gene_numbers[-1]
    last_in_sequence_gene_number = int(last_in_sequence[1:])
    print("last in sequence", last_in_sequence)

    ending_h8 = parse_grn_float2str(ending_h8_float)
    print("ending h8 float", ending_h8_float)
    last_grn = list(aligned_grns.values())[-1]
    last_grn_float = parse_grn_str2float(last_grn)
    print("aligned_grns", aligned_grns)
    last_gene_number = list(aligned_grns.keys())[-1]
    last_gene_number_int = int(last_gene_number[1:])

    missing_grns = last_in_sequence_gene_number - last_gene_number_int
    print("missing_grns", missing_grns)

    missing_h8_ = int(100 * (ending_h8_float - last_grn_float))
    missing_h8 = min(missing_grns, missing_h8_)
    missing_ctail = max(0, missing_grns - missing_h8)
    print("missing_h8", missing_h8)
    print("missing_ctail", missing_ctail)
    if '8' not in ending_h8:
        missing_ctail += missing_h8

    c_tail_float = [(100 + i + 1) for i in range(missing_ctail)]
    h8_float = grns_float[grns_float.index(last_grn_float) + 1:grns_float.index(last_grn_float) + missing_h8 + 1]
    print(h8_float)
    c_tail_float = h8_float + c_tail_float
    c_tail_float = sorted(c_tail_float)

    c_tail_str = [parse_grn_float2str(x) for x in c_tail_float]
    c_tail_list = list(zip(missing_gene_numbers[-(query_gene_len - last_gene_number_int):], c_tail_str))

    return c_tail_list, last_gene_number_int"""