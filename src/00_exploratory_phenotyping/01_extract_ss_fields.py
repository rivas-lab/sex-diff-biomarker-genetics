import pandas as pd

traits = ['2724', '3581', '3591', '2834', '3700', '3710', '3720', '2804', '3546']
f_trait = ["f.%s.%s.0" %(trait, i) for i in range(0,3) for trait in traits]
fields = ["f.eid"]+f_trait
my_f = "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/2000269/21730/download/ukb21730.tab"
tab = pd.read_csv(my_f, sep="\t", usecols=fields)
tab.to_csv("../../data/hormone_med/01_extract/ss_dat.tsv", sep="\t", index=False)

preg_f = "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/21732/download/ukb21732.tab"
preg_list = ["f.3140.%s.0" %(i) for i in range(0,3)]
preg_fields = ["f.eid"]+preg_list
preg_tab = pd.read_csv(preg_f, sep="\t", usecols=preg_fields)
preg_tab.to_csv("../../data/hormone_med/01_extract/preg_dat.tsv", sep="\t", index=False)
