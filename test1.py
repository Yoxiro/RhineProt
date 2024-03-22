import RhineAMP
Protein_Sequence = "LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTESA"
dpc = RhineAMP.DPC(Protein_Sequence)
cksaap = RhineAMP.CKSAAP(Protein_Sequence,k=0)

for i in dpc:
    print(i)