import pandas as pd
ensg = pd.read_csv("./mm10_ensembl_v91_ensg.tsv", sep = "\t", header = None, names = ["A", "B", "C", "D", "E"], index_col = 3)

data = pd.read_csv("./mouse.csv", index_col = 0)
data_ensg = data.merge(ensg[["D"]], right_index=True, left_index=True)
data_ensg.drop_duplicates("D", inplace=True)

data_ensg.index = data_ensg["D"]

data_ensg.drop("D", axis=1, inplace=True)
data_ensg.to_csv("./mouse_symbol.csv", index_label=False)