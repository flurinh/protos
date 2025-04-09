from Bio.Phylo.BaseTree import Tree, Clade
from Bio import Phylo
import pandas as pd


def gen_tree_rec(fam_df: pd.DataFrame, leaf_df: pd.DataFrame, idf: tuple = (), limit_depth: int = 4, max_depth : int = 4):
    level = len(idf)
    length = level
    if level > 4:
        pass
    else:
        if level == 0:
            # ROOT
            # print("initializing root")
            l1_children = list(set(fam_df['L1'].to_list()))
            clades = [gen_tree_rec(fam_df, leaf_df, (i,), limit_depth=limit_depth, max_depth=max_depth) for i in l1_children]
            return Clade(branch_length=length, name='root', clades=clades)
        elif level < max_depth:
            # NODES
            idf_ = [*idf]
            fill_ = [0 for _ in range((max_depth-1) - len(idf_))]
            idf_ += fill_
            name = fam_df[(fam_df['L1']==idf_[0]) &
                          (fam_df['L2']==idf_[1]) &
                          (fam_df['L3']==idf_[2])].reset_index().loc[0, 'class']
            # find all children of idf:
            # CASE 1: we have another column with identifiers
            if level < (max_depth - 1):
                fam_cols = fam_df.columns.to_list()
                filter_cols = fam_cols[:level]
                new_col = fam_cols[level]
                # print("Level {}:  |  family columns: {}  |  filter columns: {}  |  column of subfamily: {}".format(
                #    level, fam_cols, filter_cols, new_col))
                df_ = fam_df
                for c, col in enumerate(filter_cols):
                    df_ = df_[df_[col] == idf_[c]]
                # remove first row
                children = list(set(list(df_.loc[1:, new_col])))
            # CASE 2: we need to add children from the leaf table
            else:
                leaf_children_df =  leaf_df[(leaf_df['family_id_0']==idf[0]) &
                                            (leaf_df['family_id_1']==idf[1]) &
                                            (leaf_df['family_id_2']==idf[2])]
                children = [i+1 for i in range(len(leaf_children_df)-1)]
            # return leaf if we reached limit_depth => if limit is 2, we want to return an empty node reaching this depth
            # print("checking for depth limitations: level = {}  |  limit = {}".format(level, limit_depth))
            if limit_depth <= level:
                # print("reached depth limit, returning current node as leaf")
                return Clade(branch_length=length, name=name, width=0.1)
            # return another level
            else:
                # print("not yet", (*idf, 'X'), children, name)
                clades = [gen_tree_rec(fam_df, leaf_df, idf=(*idf, child), limit_depth=limit_depth, max_depth=max_depth)\
                          for child in children]
                return Clade(branch_length=length, name=name, clades=clades, width=.1)
        else:
            # LEAFES
            leaf = leaf_df[(leaf_df['family_id_0']==idf[0]) &
                           (leaf_df['family_id_1']==idf[1]) &
                           (leaf_df['family_id_2']==idf[2])]
            name=leaf['gene'].astype(str).iloc[idf[3]]
            return Clade(branch_length=length, name=name, width=0.1)


if __name__ == '__main__':
    fam_df_path = 'data/gpcrdb/fam_df.csv'
    seq_table_path = 'data/seq_table.csv'
    fam_tree_path = 'data/gpcrdb/fam_tree.tree'
    fam_df = pd.read_csv(fam_df_path)
    seq_table = pd.read_csv(seq_table_path)
    tree = gen_tree_rec(fam_df=fam_df, leaf_df=seq_table, idf=(), limit_depth=4)
    with open(fam_tree_path, 'w') as f:
        Phylo.NewickIO.write(tree, f)
