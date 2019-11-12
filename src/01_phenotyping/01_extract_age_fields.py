import pandas as pd
fields = ["f.eid", "f.53.0.0", "f.53.1.0", "f.53.2.0", "f.34.0.0", "f.21022.0.0", "f.21003.0.0", "f.21003.1.0", "f.21003.2.0"]
my_f = "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/24611/download/ukb24611.tab"
tab = pd.read_csv(my_f, sep="\t", usecols=fields)
tab.to_csv("../../data/hormone_med/01_extract/ukb24611_age.tsv", sep="\t", index=False)

# 2
my_f2 = "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/9797/download/ukb9797.tab"
tab2 = pd.read_csv(my_f2, sep="\t", usecols=fields)
tab2.to_csv("../../data/hormone_med/01_extract/ukb9797_age.tsv", sep="\t", index=False)

# 3
my_f3 = "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/21732/download/ukb21732.tab"
tab3 = pd.read_csv(my_f3, sep="\t", usecols=fields)
tab3.to_csv("../../data/hormone_med/01_extract/ukb21732_age.tsv", sep="\t", index=False)
