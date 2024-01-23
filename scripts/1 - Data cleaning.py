import pandas as pd
import numpy as np

# load csv file into pandas dataframe
df = pd.read_csv('CKO_bulk.csv', index_col=0)

# remove entries marked as "_" or "-" in the first column
df = df[df.index != "_"]
df = df[df.index != "-"]

# drop the duplicated gene names and keep the one with the highest read count
df = df.groupby(df.index).max()

# keep all the genes with reads > 1
# This is done to help normalizing the data with Log(10)
df = df[df.iloc[:, 0] > 1]

# convert read counts to Log(10)
df.iloc[:, 0] = np.log10(df.iloc[:, 0])

# save cleaned data to a new csv file
df.to_csv('Clean_CKO_bulk.csv')
