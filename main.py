import RhineAMP

if __name__ == "__main__":
    sequences = RhineAMP.Sequences()
    sequences.LoadFromFasta("AMPdata/uniprotkb_clean_cdhit_remove.fasta")
    sequences.AAC(core=4)

    # fa = RhineAMP.FastaRead("AMPdata/uniprot_clean_cdhit.fasta")
    # l = fa.shape[0]
    # lis = []
    # for i in range(l):
    #     lis.append(RhineAMP.FEGS(fa.iloc[i,1],header=fa.iloc[i,0]))
    # df = pandas.concat(lis)
    # df.to_csv("1.csv")
