import RhineAMP

"""
import os 
os.chdir(r"F:\code\RhineAMP")
"""

fa=FastaRead(Filepath="AMPdata/general_amps.fasta")

fa.iloc[:,2]>10

fa = fa[fa.iloc[:,2]>10]
fa = fa[fa.iloc[:,2]<100]
fa = fa[~fa['sequence'].str.contains("B")]
fa = fa[~fa['sequence'].str.contains("J")]
fa = fa[~fa['sequence'].str.contains("O")]
fa = fa[~fa['sequence'].str.contains("U")]
fa = fa[~fa['sequence'].str.contains("X")]
fa = fa[~fa['sequence'].str.contains("Z")]
