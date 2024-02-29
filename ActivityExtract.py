import os
os.chdir(r"F:\code\RhineAMP")

import pandas
fileName = r"AMPdata\general_amps.txt"
tsv_file = pandas.read_csv(
    fileName,
    sep='\t',
    header=0,
    encoding="utf-8"
)

activity = tsv_file[["DRAMP_ID","Sequence","Activity"]]

activity_type = set()
count = 0
for i in range(activity.shape[0]):
    if "Antibacterial" not in activity.iloc[i,2]:
        count +=1

