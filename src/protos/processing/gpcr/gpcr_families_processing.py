import os
import pandas as pd


def get_indent(line_text: str):
    idx = 0
    while line_text[idx] == ' ':
        idx+=1
    return idx


def get_level(indent=0):
    assert indent % 4 == 0, print("Found invalid indent")
    return indent // 4


def parse_line(line_text: str):
    options = line_text.split(',\"')
    if len(options) > 10:
        #print(options)
        gene = options[10].replace('\"', '')
        if len(gene) == 0:
            gene = options[9]
            print(options)
        fullname = options[11].replace('\"', '')
        uniprot = options[15].replace('\"', '')
        return [fullname, gene, uniprot]
    else:
        return line_text.replace.replace('\n', '')


def get_fam_identifier(prev_fam_identifier: list[int] = [-1, -1, -1, -1], new_index: int = 0):
    fam_identifier = prev_fam_identifier
    reset = False
    for j in range(len(prev_fam_identifier)):
        if reset:
            fam_identifier[j] = 0
        elif j == new_index:
            # we are in the class of the new index
            # go to next class
            fam_identifier[j] = prev_fam_identifier[j] + 1
            # reset all sub-class indices
            if j < len(prev_fam_identifier)-1:
                reset=True
        else:
            continue
    return fam_identifier


def map_fam_ids(fam_df, id0, id1, id2):
    f1 = fam_df[(fam_df['L1']==id0) & (fam_df['L2']==0) & (fam_df['L3']==0)].iloc[0]['class']
    f2 = fam_df[(fam_df['L1']==id0) & (fam_df['L2']==id1) & (fam_df['L3']==0)].iloc[0]['class']
    f3 = fam_df[(fam_df['L1']==id0) & (fam_df['L2']==id1) & (fam_df['L3']==id2)].iloc[0]['class']
    return f1, f2, f3


if __name__ == '__main__':
    families_path = 'data/gpcrdb/families/families.txt'
    fam_df = pd.DataFrame([], columns = ['L1', 'L2', 'L3', 'class']) # dict
    with open(families_path) as f:
        family_dict = {}
        gene_data = []
        lines = f.readlines()
        prev_fam_identifier = [0, -1, -1, -1]
        prev_indent = 0
        for l, line in enumerate(lines):
            indent = get_indent(line)
            level = get_level(indent)
            fam_identifier = get_fam_identifier(prev_fam_identifier, level)
            if level <= 2:
                fam_df.loc[len(fam_df)] = fam_identifier[:-1] + [line[indent:]]
            if level == 3:
                # add data of gene
                gene_info = parse_line(line[indent:])
                gene_data.append(gene_info + fam_identifier[:-1])
            prev_fam_identifier = fam_identifier

    df = pd.DataFrame(gene_data, columns=['name', 'gene', 'uniprot', 'family_id_0', 'family_id_1', 'family_id_2'])

    df[['f1', 'f2', 'f3']] = df.apply(lambda x: map_fam_ids(fam_df, x.family_id_0, x.family_id_1, x.family_id_2), axis=1,
                                      result_type='expand')

    fam_df.to_csv('data/gpcrdb/fam_df.csv', index=False)
    df.to_csv('data/seq_table.csv', index=False)