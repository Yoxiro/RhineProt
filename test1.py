import RhineAMP

fa = RhineAMP.FastaRead("AMPdata/uniprotkb_clean_cdhit_remove.fasta")

fa = RhineAMP.FastaPadding(fa)

fa = RhineAMP.BLOSUM_62_fasta(fa)

print(fa)