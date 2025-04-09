


def get_occs_of_cys(grn_table, grn_start='8x52'):
    # count number of cysteins in c tail
    columns = grn_table.columns.tolist()
    idx = columns.index(grn_start)
    cols = columns[idx:]
    results = {}
    indices = grn_table.index.tolist()
    for i, idx in enumerate(indices):
        occs = 0
        for g, grn in enumerate(cols):
            x = grn_table.loc[idx, grn]
            if not isinstance(x, str):
                x = x.iloc[0]
            if x != '-':
                if 'C' in x:
                    occs += 1
        results.update({idx: occs})
    return results


def get_dists_from_cterminus(grn_table, grn_start='8x52'):
    # get distance to closest cys starting from the c terminus
    columns = grn_table.columns.tolist()
    idx = columns.index(grn_start)
    cols = columns[idx:][::-1]  # all grns from end to start (reverse order)
    results = {}
    indices = grn_table.index.tolist()
    for i, idx in enumerate(indices):
        res = []
        start = 0
        start_uidx = None
        for g, grn in enumerate(cols):
            x = grn_table.loc[idx, grn]

            if not isinstance(x, str):
                x = x.iloc[0]
            if (x != '-') & (start == 0):
                start = cols.index(grn)
                start_uidx = int(x[1:])
                if 'C' in x:
                    res.append(0)
            elif x != '-':
                if 'C' in x:
                    res.append(start_uidx-int(x[1:]))
        results.update({idx: res})
    return results


def get_dists_from_h8(grn_table, grn_start='8x52'):
    # get distance to first cys starting from the h8 helix
    columns = grn_table.columns.tolist()
    idx = columns.index(grn_start)
    cols = columns[idx:]
    results = {}
    indices = grn_table.index.tolist()
    for i, idx in enumerate(indices):
        res = []
        start = 0
        for g, grn in enumerate(cols):
            x = grn_table.loc[idx, grn]
            if x != '-':
                if g == 0:
                    start = int(x[1:])
                if 'C' in x:
                    res.append(int(x[1:])-start)
        if len(res) >= 1:
            results.update({idx: res})
        else:
            results.update({idx: [0]})
    return results