# RhineProt
Under development

# Usage

Download the code, and makesure the **RhineAMP** folder is under the your project root

## FEGS module

To use the FEGS module:

with one Protein Sequence:
```
import RhineAMP

fegs = RhineAMP.FEGS(Protein_Sequence,header)
```

fegs contains 580 features in pandas.Series form and its name is the corresponding header of Sequence

If there are many protein sequences that need to extract feature values, you can use multiprocessing computation to save time:
```
import RhineAMP

sequences = RhineAMP.Sequences()
sequences.LoadFromFasta(filepath)
sequences.FEGS(core = 4, save_result= True, filepath="Result.csv")
```