import pandas


def FastaFilter(Fasta: pandas.DataFrame,
                Filter_words: list = None) -> pandas.DataFrame:
    if Filter_words is None:
        Filter_words = ["anti", "Anti"]
    for char in Filter_words:
        Fasta = Fasta[~Fasta["header"].str.contains(char)]
    return Fasta
