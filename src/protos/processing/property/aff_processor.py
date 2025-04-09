import pandas as pd
import numpy as np


class AffinityProcessor:
    def __init__(self,
                 folder='data/gpcrdb/affinities/',
                 filename='affinities_B.csv',
                 load_data=True,
                 process_xls=False,
                 xls_path='data/gpcrdb/couplings/gpcr_gprotein_couplings.xlsx'):
        self.folder = folder
        self.xls_path = xls_path
        self.sheets = None
        self.xls = None
        self.data = pd.DataFrame()
        if process_xls:
            self.xls = pd.ExcelFile(self.xls_path)
            self.sheets = self.xls.sheet_names
        if load_data:
            self.load_data(filename)

    def sheet2csv(self, sheetname='B-import', filename='affinities_B.csv'):
        if self.sheets == None:
            self.xls = pd.ExcelFile(self.xls_path)
            self.sheets = self.xls.sheet_names
        if sheetname not in self.sheets:
            print("Invalid sheetname, using default: ", self.sheets[0])
            sheetname = self.sheets[0]
        df = pd.read_excel(self.path, sheet_name=sheetname)
        columns = df.iloc[0]
        columns = [' '.join(col.split('\n')) for col in columns]
        df.columns = columns
        df['uni_entry_name'] = df['Receptor (UniProt)'].apply(lambda x: x + '_HUMAN')
        df.to_csv(self.folder+filename, index=False)

    def load_data(self, filename='affinities_B.csv'):
        self.data = pd.read_csv(self.folder+filename)

    def get_prim_prot(self):
        cols = ['log(Emax /EC50) Gs', 'log(Emax /EC50) Gi/o', 'log(Emax /EC50) Gq/11', 'log(Emax /EC50) G12/13']
        self.data['pref'] = self.data[cols].idxmax(axis=1)
        self.data['pref'] = self.data['pref'].apply(lambda x: x.split(' ')[-1])

    def get_max_sel(self):
        cols = ['log(Emax /EC50) Gs', 'log(Emax /EC50) Gi/o', 'log(Emax /EC50) Gq/11', 'log(Emax /EC50) G12/13']
        self.data['max'] = self.data[cols].max(axis=1)

    def get_diff_sel(self):
        cols = ['log(Emax /EC50) Gs', 'log(Emax /EC50) Gi/o', 'log(Emax /EC50) Gq/11', 'log(Emax /EC50) G12/13']
        df = self.data[cols]
        df = pd.DataFrame(np.sort(df.values)[:, -2:], columns=['second', 'largest'])
        self.data['differential_prim_sec'] = df['second'] / df['largest']
        del df