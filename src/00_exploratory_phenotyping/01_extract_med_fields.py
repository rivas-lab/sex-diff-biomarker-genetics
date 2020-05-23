import pandas as pd

med_cols = ["f.20003.%d.%d" %(i,j) for i in range(0,3) for j in range(0,48)]
fields = ["f.eid"] + med_cols
my_f = "/oak/stanford/groups/mrivas/ukbb/24983/phenotypedata/9796/24611/download/ukb24611.tab"
tab = pd.read_csv(my_f, sep="\t", usecols=fields)
tab.to_csv("../../data/hormone_med/01_extract/med_extract.tsv", sep="\t", index=False)
