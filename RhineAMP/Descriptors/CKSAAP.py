import sys

import pandas

def CKSAAP(Protein_Sequence:str,k:int,header:str = None)->pandas.Series:
    """

    :param Protein_Sequence:
    :param k: k is the distance between residue pairs
    :param header: the header of the given protein sequence(default:"0")
    :return:
    """
    """
    Protein_Sequence = "ABCDEFG"
    Protein_Length = 7
    k = 2 
    then cksaap.len = 4 max_index = 3 = 7 - 1 - 2
    k = 1
    then cksaap.len = 5 max_index = 4 = 7 - 1 - 2
    k = 0
    then cksaap.len = 6
    """
    if header is None:
        header = "0"
    Protein_Sequence = Protein_Sequence.strip()
    Protein_Length = len(Protein_Sequence)
    if k >= Protein_Length:
        print("the length of the Protein should be greater than the k")
        sys.exit(0)

    index_string = ["A", "I", "L", "V", "M", "F", "W", "Y",
                    "N", "C", "Q", "S", "T",
                    "R", "H", "K",
                    "D", "E",
                    "G", "P"]

    index_string_after = [i + j for i in index_string for j in index_string]
    cksaap = pandas.Series([0 for _ in index_string_after],
                               index=index_string_after,
                               dtype=float,
                               name=header)
    for char_index in range(Protein_Length - k - 1):
        chars = Protein_Sequence[char_index:char_index + 2]
        if chars not in index_string_after:
            raise Exception("Invalid amino acid symbols")
        cksaap[chars] += 1
    cksaap = cksaap / (Protein_Length - k - 1)
    return cksaap