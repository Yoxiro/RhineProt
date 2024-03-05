import pandas


def FastaClean(Fasta: pandas.DataFrame,
               LengthLimit: tuple or list = (10, 100),
               Chars=None):
    if Chars is None:
        Chars = ["B", "J", "O", "U", "X", "Z"]
    minimum = LengthLimit[0]
    maximum = LengthLimit[1]
    Fasta = Fasta[Fasta.iloc[:, 2] >= minimum]
    Fasta = Fasta[Fasta.iloc[:, 2] <= maximum]
    for char in Chars:
        Fasta = Fasta[~Fasta["sequence"].str.contains(char)]
    return Fasta
