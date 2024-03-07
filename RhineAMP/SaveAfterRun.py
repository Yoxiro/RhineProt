import pandas
def SaveAfterRun(Func,filepath:str,Protein_Sequence):
    res = Func(Protein_Sequence)
    res.to_csv(filepath)